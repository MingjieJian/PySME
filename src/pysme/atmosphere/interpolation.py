# -*- coding: utf-8 -*-
""" Handles reading and interpolation of atmopshere (grid) data """
import itertools
import logging

import numpy as np
from astropy import constants as const
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from ..large_file_storage import setup_atmo
from .atmosphere import Atmosphere as Atmo
from .atmosphere import AtmosphereError, AtmosphereGrid
from .savfile import SavFile

logger = logging.getLogger(__name__)

# Radius and Surface Gravity of the Sun
# Required for spherical models
R_sun = const.R_sun.to_value("cm")
g_sun = (const.G * const.M_sun / const.R_sun ** 2).to_value("cm/s**2")
logg_sun = np.log10(g_sun)


class AtmosphereInterpolator:
    def __init__(self, depth=None, interp=None, geom=None, lfs_atmo=None, verbose=0):
        self.depth = depth
        self.interp = interp
        self.geom = geom
        if lfs_atmo is None:
            lfs_atmo = setup_atmo()
        self.lfs_atmo = lfs_atmo

        self.source = None
        self.atmo_grid = None
        self.verbose = verbose

    def interp_atmo_grid(self, atmo_grid, teff, logg, monh):
        """
        General routine to interpolate in 3D grid of model atmospheres

        Parameters
        -----
        Teff : float
            effective temperature of desired model (K).
        logg : float
            logarithmic gravity of desired model (log cm/s/s).
        MonH : float
            metalicity of desired model.
        atmo_in : Atmosphere
            Input atmosphere
        verbose : {0, 1}, optional
            how much information to plot/print (default: 0)
        plot : bool, optional
            wether to plot debug information (default: False)
        reload : bool
            wether to reload atmosphere information from disk (default: False)

        Returns
        -------
        atmo : Atmosphere
            interpolated atmosphere data
        """

        # Internal parameters.
        if not isinstance(atmo_grid, AtmosphereGrid):
            if self.atmo_grid is None or self.source != atmo_grid:
                atmo_file = self.lfs_atmo.get(atmo_grid)
                self.source = atmo_grid
                self.atmo_grid = SavFile(
                    atmo_file, source=self.source, lfs=self.lfs_atmo
                )
        else:
            self.atmo_grid = atmo_grid
            self.source = atmo_grid.source

        if self.geom is not None:
            if self.geom == "PP":
                atmo_grid = self.atmo_grid[self.atmo_grid.radius <= 1]
            elif self.geom == "SPH":
                atmo_grid = self.atmo_grid[self.atmo_grid.radius > 1]
        else:
            atmo_grid = self.atmo_grid

        assert len(atmo_grid) > 0, "No atmospheres for this geometry"

        # Get field names in ATMO and ATMO_GRID structures.
        depth = self.determine_depth_scale(self.depth, atmo_grid)
        interp = self.determine_interpolation_scale(self.interp, atmo_grid)

        # Find the corner models bracketing the given values
        icor = self.find_corner_models(teff, logg, monh, atmo_grid)

        # Interpolate the corner models
        atmo = self.interpolate_corner_models(
            teff, logg, monh, icor, atmo_grid, interp=interp
        )

        # TODO: Or should we only consider spherical models for interpolation of spherical modesl are requested?
        geom, radius = self.spherical_model_correction(atmo_grid, icor, logg)
        # Create ATMO.GEOM, if necessary, and set value.
        if self.geom is not None and self.geom != geom:
            if self.geom == "SPH":
                raise AtmosphereError(
                    "Input ATMO.GEOM='%s' was requested but the model only supports PP (at this point)."
                    % self.geom
                )
            else:
                logger.info(
                    "Input ATMO.GEOM='%s' overrides '%s' from grid.",
                    self.geom,
                    geom,
                )

        # Add standard ATMO input fields, if they are missing from ATMO_IN.
        atmo.depth = depth
        atmo.interp = interp
        atmo.geom = geom
        if radius is not None:
            atmo.radius = radius
        atmo.source = self.source
        atmo.method = "grid"

        return atmo

    def interp_atmo_pair(self, atmo1, atmo2, frac, interpvar="RHOX", itop=0):
        """
        Interpolate between two model atmospheres, accounting for shifts in
        the mass column density or optical depth scale.

        How it works:

        1. The second atmosphere is fitted onto the first, individually for
        each of the four atmospheric quantitites: T, xne, xna, rho.
        The fitting uses a linear shift in both the (log) depth parameter and
        in the (log) atmospheric quantity. For T, the midpoint of the two
        atmospheres is aligned for the initial guess. The result of this fit
        is used as initial guess for the other quantities. A penalty function
        is applied to each fit, to avoid excessively large shifts on the
        depth scale.
        2. The mean of the horizontal shift in each parameter is used to
        construct the common output depth scale.
        3. Each atmospheric quantity is interpolated after shifting the two
        corner models by the amount determined in step 1), rescaled by the
        interpolation fraction (frac).

        Parameters
        ------
        atmo1 : Atmosphere
            first atmosphere to interpolate
        atmo2 : Atmosphere
            second atmosphere to interpolate
        frac : float
            interpolation fraction: 0.0 -> atmo1 and 1.0 -> atmo2
        interpvar : {"RHOX", "TAU"}, optional
            atmosphere interpolation variable (default:"RHOX").
        itop : int, optional
            index of top point in the atmosphere to use. default
            is to use all points (itop=0). use itop=1 to clip top depth point.
        atmop : array[5, ndep], optional
            interpolated atmosphere prediction (for plots)
            Not needed if atmospheres are provided as structures. (default: None)
        verbose : {0, 5}, optional
            diagnostic print level (default 0: no printing)
        plot : {-1, 0, 1}, optional
            diagnostic plot level. Larger absolute
            value yields more plots. Negative values cause a wait for keypress
            after each plot. (default: 0, no plots)
        old : bool, optional
            also plot result from the old interpkrz2 algorithm. (default: False)

        Returns
        ------
        atmo : Atmosphere
            interpolated atmosphere
            .RHOX (vector[ndep]) mass column density (g/cm^2)
            .TAU  (vector[ndep]) reference optical depth (at 5000 Å)
            .temp (vector[ndep]) temperature (K)
            .xne  (vector[ndep]) electron number density (1/cm^3)
            .xna  (vector[ndep]) atomic number density (1/cm^3)
            .rho  (vector[ndep]) mass density (g/cm^3)
        """
        # """
        # History
        # -------
        # 2004-Apr-15 Valenti
        #     Initial coding.
        # MB
        #     interpolation on TAU scale
        # 2012-Oct-30 TN
        #     Rewritten to use either column mass (RHOX) or
        #     reference optical depth (TAU) as vertical scale. Shift-interpolation
        #     algorithms have been improved for stability in cool dwarfs (<=3500 K).
        #     The reference optical depth scale is preferred in terms of interpolation
        #     accuracy across most of parameter space, with significant improvement for
        #     both cool models (where depth vs temperature is rather flat) and hot
        #     models (where depth vs temperature exhibits steep transitions).
        #     Column mass depth is used by default for backward compatibility.
        # 2013-May-17 Valenti
        #     Use frac to weight the two shifted depth scales,
        #     rather than simply averaging them. This change fixes discontinuities
        #     when crossing grid nodes.
        # 2013-Sep-10 Valenti
        #     Now returns an atmosphere structure instead of a
        #     [5,NDEP] atmosphere array. This was necessary to support interpolation
        #     using one variable (e.g., TAU) and radiative transfer using a different
        #     variable (e.g. RHOX). The atmosphere array could only store one depth
        #     variable, meaning the radiative transfer variable had to be the same
        #     as the interpolation variable. Returns atmo.RHOX if available and also
        #     atmo.TAU if available. Since both depth variables are returned, if
        #     available, this routine no longer needs to know which depth variable
        #     will be used for radiative transfer. Only the interpolation variable
        #     is important. Thus, the interpvar= keyword argument replaces the
        #     type= keyword argument. Very similar code blocks for each atmospheric
        #     quantity have been unified into a single code block inside a loop over
        #     atmospheric quantities.
        # 2013-Sep-21 Valenti
        #     Fixed an indexing bug that affected the output depth
        #     scale but not other atmosphere vectors. Itop clipping was not being
        #     applied to the depth scale ('RHOX' or 'TAU'). Bug fixed by adding
        #     interpvar to vtags. Now atmospheres interpolated with interp_atmo_grid
        #     match output from revision 398. Revisions back to 399 were development
        #     only, so no users should be affected.
        # 2014-Mar-05 Piskunov
        #     Replicated the removal of the bad top layers in models
        #     for interpvar eq 'TAU'
        # """

        # Internal program parameters.
        min_drhox = min_dtau = 0.01  # minimum fractional step in RHOX
        # min_dtau = 0.01  # minimum fractional step in TAU
        interpvar = interpvar.lower()
        ##
        ## Select interpolation variable (RHOX vs. TAU)
        ##

        # Check which depth scales are available in both input atmospheres.
        tags1 = atmo1.dtype.names
        tags2 = atmo2.dtype.names
        ok_tau = "tau" in tags1 and "tau" in tags2
        ok_rhox = "rhox" in tags1 and "rhox" in tags2
        if not ok_tau and not ok_rhox:
            raise AtmosphereError(
                "atmo1 and atmo2 structures must both contain RHOX or TAU"
            )

        # Set interpolation variable, if not specified by keyword argument.
        if interpvar is None:
            interpvar = "tau" if ok_tau else "rhox"
        if interpvar != "tau" and interpvar != "rhox":
            raise AtmosphereError("interpvar must be 'TAU' (default) or 'RHOX'")

        has_spherical_height = (
            getattr(atmo1, "radius", 0) > 1
            and getattr(atmo2, "radius", 0) > 1
            and "height" in tags1
            and "height" in tags2
        )

        def to_interp_space(vtag, values):
            values = np.asarray(values)
            if vtag == "height":
                return values
            return np.log10(values)

        def from_interp_space(vtag, values):
            values = np.asarray(values)
            if vtag == "height":
                return values
            return 10.0 ** values

        ##
        ## Define depth scale for both atmospheres
        ##

        # Define depth scale for atmosphere #1
        mask1 = np.full(len(atmo1[interpvar]), True)
        mask1[:itop] = False

        itop1 = itop
        while (atmo1[interpvar][itop1 + 1] / atmo1[interpvar][itop1] - 1) <= min_drhox:
            mask1[itop1] = False
            itop1 += 1

        mask1 &= atmo1[interpvar] != 0
        mask1 &= np.isfinite(atmo1[interpvar])
        ndep1 = np.count_nonzero(mask1)

        depth1 = np.log10(atmo1[interpvar][mask1])

        # Define depth scale for atmosphere #2
        mask2 = np.full(len(atmo2[interpvar]), True)
        mask2[:itop] = False
        itop2 = itop
        while (atmo2[interpvar][itop2 + 1] / atmo2[interpvar][itop2] - 1) <= min_drhox:
            mask2[itop2] = False
            itop2 += 1

        mask2 &= atmo2[interpvar] != 0
        mask2 &= np.isfinite(atmo2[interpvar])
        ndep2 = np.count_nonzero(mask2)

        depth2 = np.log10(atmo2[interpvar][mask2])

        ##
        ## Prepare to find best shift parameters for each atmosphere vector.
        ##

        # List of atmosphere vectors that need to be shifted.
        # The code below assumes 'TEMP' is the first vtag in the list.
        vtags = ["temp", "xne", "xna", "rho", interpvar]
        if interpvar == "rhox" and ok_tau:
            vtags += ["tau"]
        if interpvar == "tau" and ok_rhox:
            vtags += ["rhox"]
        if has_spherical_height:
            vtags += ["height"]
        nvtag = len(vtags)

        # Adopt arbitrary uncertainties for shift determinations.
        err1 = np.full(ndep1, 0.05)

        # Initial guess for TEMP shift parameters.
        # Put depth and TEMP midpoints for atmo1 and atmo2 on top of one another.
        npar = 4
        ipar = np.zeros(npar, dtype="f4")
        temp1 = np.log10(atmo1.temp[mask1])
        temp2 = np.log10(atmo2.temp[mask2])
        mid1 = np.argmin(np.abs(temp1 - 0.5 * (temp1[1] + temp1[-2])))
        mid2 = np.argmin(np.abs(temp2 - 0.5 * (temp2[1] + temp2[-2])))
        ipar[0] = depth1[mid1] - depth2[mid2]  # horizontal
        ipar[1] = temp1[mid1] - temp2[mid2]  # vertical

        # Apply a constraint on the fit, to avoid runaway for cool models, where
        # the temperature structure is nearly linear with both TAU and RHOX.
        constraints = np.zeros(npar)
        constraints[0] = 0.5  # weakly constrain the horizontal shift

        # For first pass ('TEMP'), use all available depth points.
        igd = np.isfinite(depth1)
        ngd = igd.size

        ##
        ## Find best shift parameters for each atmosphere vector.
        ##

        # Loop through atmosphere vectors.
        pars = np.zeros((nvtag, npar))
        for ivtag, vtag in enumerate(vtags):

            # Find vector in each structure.
            if vtag not in tags1:
                raise AtmosphereError("atmo1 does not contain " + vtag)
            if vtag not in tags2:
                raise AtmosphereError("atmo2 does not contain " + vtag)

            vect1 = to_interp_space(vtag, atmo1[vtag][mask1])
            vect2 = to_interp_space(vtag, atmo2[vtag][mask2])

            # Fit the second atmosphere onto the first by finding the best horizontal
            # shift in depth2 and the best vertical shift in vect2.
            pars[ivtag], _ = self.interp_atmo_constrained(
                depth1[igd],
                vect1[igd],
                err1[igd],
                ipar,
                x2=depth2,
                y2=vect2,
                y1=vect1,
                ndep=ngd,
                constraints=constraints,
            )

            # After first pass ('TEMP'), adjust initial guess and restrict depth points.
            if ivtag == 0:
                ipar = [pars[0, 0], 0.0, 0.0, 0.0]
                try:
                    igd = np.where(
                        (depth1 >= min(depth2[igd]) + ipar[0])
                        & (depth1 <= max(depth2[igd]) + ipar[0])
                    )[0]
                except IndexError:
                    pass
                if igd.size < 2:
                    raise AtmosphereError("unstable shift in temperature")

        ##
        ## Use mean shift to construct output depth scale.
        ##

        # Calculate the mean depth2 shift for all atmosphere vectors.
        xsh = np.sum(pars[:, 0]) / nvtag

        # Base the output depth scale on the input scale with the fewest depth points.
        # Combine the two input scales, if they have the same number of depth points.
        depth1f = depth1 - xsh * frac
        depth2f = depth2 + xsh * (1 - frac)
        if ndep1 > ndep2:
            depth = depth2f
        elif ndep1 == ndep2:
            depth = depth1f * (1 - frac) + depth2f * frac
        elif ndep1 < ndep2:
            depth = depth1f
        ndep = len(depth)

        ##
        ## Interpolate input atmosphere vectors onto output depth scale.
        ##

        # Loop through atmosphere vectors.
        vects = np.zeros((nvtag, ndep))
        for ivtag, (vtag, par) in enumerate(zip(vtags, pars)):

            # Extract data
            vect1 = to_interp_space(vtag, atmo1[vtag][mask1])
            vect2 = to_interp_space(vtag, atmo2[vtag][mask2])

            # Identify output depth points that require extrapolation of atmosphere vector.
            depth1f = depth1 - par[0] * frac
            depth2f = depth2 + par[0] * (1 - frac)
            x1max = np.max(depth1f)
            x2max = np.max(depth2f)
            iup = (depth > x1max) | (depth > x2max)
            nup = np.count_nonzero(iup)
            checkup = (nup >= 1) and abs(frac - 0.5) <= 0.5 and ndep1 == ndep2

            # Combine shifted vect1 and vect2 structures to get output vect.
            vect1f = self.interp_atmo_func(depth, -frac * par, x2=depth1, y2=vect1)
            vect2f = self.interp_atmo_func(depth, (1 - frac) * par, x2=depth2, y2=vect2)
            vect = (1 - frac) * vect1f + frac * vect2f
            ends = [vect1[ndep1 - 1], vect[ndep - 1], vect2[ndep2 - 1]]
            if (
                checkup
                and np.median(ends) != vect[ndep - 1]
                and (
                    abs(vect1[ndep1 - 1] - 4.2) < 0.1
                    or abs(vect2[ndep2 - 1] - 4.2) < 0.1
                )
            ):
                vect[iup] = vect2f[iup] if x1max < x2max else vect1f[iup]
            vects[ivtag] = vect

        ##
        ## Construct output structure
        ##

        # Construct output structure with interpolated atmosphere.
        # Might be wise to interpolate abundances, in case those ever change.
        atmo = Atmo(interp=interpvar)
        stags = ["teff", "logg", "monh", "vturb", "lonh", "abund"]
        ndep_orig = len(atmo1.temp)
        for tag in tags1:

            # Default is to copy value from atmo1. Trim vectors.
            value = atmo1[tag]
            if np.size(value) == ndep_orig and tag != "abund":
                value = value[:ndep]

            # Vector quantities that have already been interpolated.
            if tag in vtags:
                ivtag = [i for i in range(nvtag) if tag == vtags[i]][0]
                value = from_interp_space(tag, vects[ivtag])

            # Scalar quantities that should be interpolated using frac.
            if tag in stags:
                if tag in tags2:
                    value = (1 - frac) * atmo1[tag] + frac * atmo2[tag]
                else:
                    value = atmo1[tag]

            # Remaining cases.
            if tag == "ndep":
                value = ndep

            # Abundances
            if tag == "abund":
                value = (1 - frac) * atmo1[tag] + frac * atmo2[tag]

            # Create or add to output structure.
            atmo[tag] = value
        return atmo

    def determine_depth_scale(self, depth, atmo_grid):
        """
        Determine ATMO.DEPTH radiative transfer depth variable. Order of precedence:

        1. Value of ATMO_IN.DEPTH, if it exists and is set
        2. Value of ATMO_GRID[0].DEPTH, if it exists and is set
        3. 'RHOX', if ATMO_GRID.RHOX exists (preferred over 'TAU' for depth)
        4. 'TAU', if ATMO_GRID.TAU exists

        Check that INTERP is valid and the corresponding field exists in ATMO.

        Parameters
        ----------
        depth : {"RHOX", "TAU", None}
            requested value, or None for autoselection based in available grid
        atmo_grid : AtmosphereGrid
            input atmosphere grid to interpolate on

        Returns
        -------
        depth : {"TAU", "RHOX"}
            The chosen depth scale

        Raises
        ------
        AtmosphereError
            If an invalid value was set in atmo_in/atmo_grid
        """

        gtags = atmo_grid.dtype.names

        if depth is not None:
            depth = depth
        elif "depth" in gtags and atmo_grid.depth is not None:
            depth = atmo_grid.depth
        elif "rhox" in gtags:
            depth = "RHOX"
        elif "tau" in gtags:
            depth = "TAU"
        else:
            raise AtmosphereError("no value for ATMO.DEPTH")
        if depth != "TAU" and depth != "RHOX":
            raise AtmosphereError(
                "ATMO.DEPTH must be 'TAU' or 'RHOX', not '%s'" % depth
            )
        if depth.lower() not in gtags:
            raise AtmosphereError(
                "ATMO.DEPTH='{}', but ATMO. {} does not exist".format(depth, depth)
            )
        return depth

    def determine_interpolation_scale(self, interp, atmo_grid):
        """
        Determine ATMO.INTERP interpolation variable. Order of precedence:

        1. Value of ATMO_IN.INTERP, if it exists and is set
        2. Value of ATMO_GRID[0].INTERP, if it exists and is set
        3. 'TAU', if ATMO_GRID.TAU exists (preferred over 'RHOX' for interpolation)
        4. 'RHOX', if ATMO_GRID.RHOX exists

        Check that INTERP is valid and the corresponding field exists in ATMO.

        Parameters
        ----------
        interp : {"RHOX", "TAU", None}
            requested interpolation axis, or None for autoselect
        atmo_grid : AtmosphereGrid
            Atmosphere grid for interpolation

        Returns
        -------
        interp: {"RHOX", "TAU"}
            the interpolation axis

        Raises
        ------
        AtmosphereError
            if a non valid parameter is set in atmo_in/atmo_grid
        """

        gtags = atmo_grid.dtype.names

        if interp is not None:
            interp = interp
        elif atmo_grid.interp is not None:
            interp = atmo_grid.interp
        elif "tau" in gtags:
            interp = "TAU"
        elif "rhox" in gtags:
            interp = "RHOX"
        else:
            raise AtmosphereError("no value for ATMO.INTERP")
        if interp not in ["TAU", "RHOX"]:
            raise AtmosphereError(
                "ATMO.INTERP must be 'TAU' or 'RHOX', not '%s'" % interp
            )
        if interp.lower() not in gtags:
            raise AtmosphereError(
                "ATMO.INTERP='{}', but ATMO. {} does not exist".format(interp, interp)
            )
        return interp

    def find_corner_models(self, teff, logg, monh, atmo_grid):
        """
        Find the models in the grid that bracket the given stellar parameters

        The purpose of the first major set of code blocks is to find values
        of [M/H] in the grid that bracket the requested [M/H]. Then in this
        subset of models, find values of log(g) in the subgrid that bracket
        the requested log(g). Then in this subset of models, find values of
        Teff in the subgrid that bracket the requested Teff. The net result
        is a set of 8 models in the grid that bracket the requested stellar
        parameters. Only these 8 "corner" models will be used in the
        interpolation that constitutes the remainder of the program.

        Parameters
        ----------
        teff : float
            effective temperature
        logg : float
            surface gravity
        monh : float
            metallicity
        atmo_grid : AtmosphereGrid
            atmosphere grid to search models in
        verbose : int, optional
            determines how many debugging messages to print, by default 0
        """

        nb = 2  # number of bracket points

        # *** DETERMINATION OF METALICITY BRACKET ***
        # Find unique set of [M/H] values in grid.
        mlist = np.unique(atmo_grid.monh)  # list of unique [M/H]

        # Test whether requested metalicity is in grid.
        mmin = np.min(mlist)  # range of [M/H] in grid
        mmax = np.max(mlist)
        if monh > mmax:  # true: [M/H] too large
            logger.info(
                "interp_atmo_grid: requested [M/H] (%.3f) larger than max grid value (%.3f). extrapolating.",
                monh,
                mmax,
            )
        if monh < mmin:  # true: logg too small
            raise AtmosphereError(
                "interp_atmo_grid: requested [M/H] (%.3f) smaller than min grid value (%.3f). returning."
                % (monh, mmin)
            )

        # Find closest two [M/H] values in grid that bracket requested [M/H].
        if monh <= mmax:
            mlo = np.max(mlist[mlist <= monh])
            mup = np.min(mlist[mlist >= monh])
        else:
            mup = mmax
            mlo = np.max(mlist[mlist < mup])
        mb = [mlo, mup]  # bounding [M/H] values

        # Trace diagnostics.
        if self.verbose >= 5:
            logger.info("[M/H]: %.3f < %.3f < %.3f", mlo, monh, mup)

        # *** DETERMINATION OF LOG(G) BRACKETS AT [M/H] BRACKET VALUES ***
        # Set up for loop through [M/H] bounds.
        gb = np.zeros((nb, nb))  # bounding gravities
        for i in range(nb):
            # Find unique set of gravities at boundary below [M/H] value.
            im = atmo_grid.monh == mb[i]  # models on [M/H] boundary
            glist = np.unique(atmo_grid[im].logg)  # list of unique gravities

            # Test whether requested logarithmic gravity is in grid.
            gmin = np.min(glist)  # range of gravities in grid
            gmax = np.max(glist)
            if logg > gmax:  # true: logg too large
                logger.info(
                    "interp_atmo_grid: requested log(g) (%.3f) larger than max grid value (%.3f). extrapolating.",
                    logg,
                    gmax,
                )

            if logg < gmin:  # true: logg too small
                raise AtmosphereError(
                    "interp_atmo_grid: requested log(g) (%.3f) smaller than min grid value (%.3f). returning."
                    % (logg, gmin)
                )

            # Find closest two gravities in Mlo subgrid that bracket requested gravity.
            if logg <= gmax:
                glo = np.max(glist[glist <= logg])
                gup = np.min(glist[glist >= logg])
            else:
                gup = gmax
                glo = np.max(glist[glist < gup])
            gb[i] = [glo, gup]  # store boundary values.

            # Trace diagnostics.
            if self.verbose >= 5:
                logger.info(
                    "log(g) at [M/H]=%.3f: %.3f < %.3f < %.3f",
                    mb[i],
                    glo,
                    logg,
                    gup,
                )

        # End of loop through [M/H] bracket values.
        # *** DETERMINATION OF TEFF BRACKETS AT [M/H] and LOG(G) BRACKET VALUES ***
        # Set up for loop through [M/H] and log(g) bounds.
        tb = np.zeros((nb, nb, nb))  # bounding temperatures
        for ig in range(nb):
            for im in range(nb):
                # Find unique set of gravities at boundary below [M/H] value.
                it = (atmo_grid.monh == mb[im]) & (
                    atmo_grid.logg == gb[im, ig]
                )  # models on joint boundary
                tlist = np.unique(atmo_grid[it].teff)  # list of unique temperatures

                # Test whether requested temperature is in grid.
                tmin = np.min(tlist)  # range of temperatures in grid
                tmax = np.max(tlist)
                if teff > tmax:  # true: Teff too large
                    raise AtmosphereError(
                        "interp_atmo_grid: requested Teff (%i) larger than max grid value (%i). returning."
                        % (teff, tmax)
                    )
                if teff < tmin:  # true: logg too small
                    logger.info(
                        "interp_atmo_grid: requested Teff (%i) smaller than min grid value (%i). extrapolating.",
                        teff,
                        tmin,
                    )

                # Find closest two temperatures in subgrid that bracket requested Teff.
                if teff > tmin:
                    tlo = np.max(tlist[tlist <= teff])
                    tup = np.min(tlist[tlist >= teff])
                else:
                    tlo = tmin
                    tup = np.min(tlist[tlist > tlo])
                tb[im, ig, :] = [tlo, tup]  # store boundary values.

                # Trace diagnostics.
                if self.verbose >= 5:
                    logger.info(
                        "Teff at log(g)=%.3f and [M/H]=%.3f: %i < %i < %i",
                        gb[im, ig],
                        mb[im],
                        tlo,
                        teff,
                        tup,
                    )

        # End of loop through log(g) and [M/H] bracket values.

        # Find and save atmo_grid indices for the 8 corner models.
        icor = np.zeros((nb, nb, nb), dtype=int)
        for it, ig, im in itertools.product(range(nb), repeat=3):
            iwhr = np.where(
                (atmo_grid.teff == tb[im, ig, it])
                & (atmo_grid.logg == gb[im, ig])
                & (atmo_grid.monh == mb[im])
            )[0]
            nwhr = iwhr.size
            if nwhr != 1:
                logger.info(
                    "interp_atmo_grid: %i models in grid with [M/H]=%.1f, log(g)=%.1f, and Teff=%i",
                    nwhr,
                    mb[im],
                    gb[im, ig],
                    tb[im, ig, it],
                )
            icor[im, ig, it] = iwhr[0]

        # Trace diagnostics.
        if self.verbose >= 1:
            logger.info("Teff=%i,  log(g)=%.3f,  [M/H]=%.3f:", teff, logg, monh)
            logger.info("indx  M/H  g   Teff     indx  M/H  g   Teff")
            for im in range(nb):
                for ig in range(nb):
                    i0 = icor[im, ig, 0]
                    i1 = icor[im, ig, 1]
                    logger.info(
                        i0,
                        atmo_grid[i0].monh,
                        atmo_grid[i0].logg,
                        atmo_grid[i0].teff,
                        i1,
                        atmo_grid[i1].monh,
                        atmo_grid[i1].logg,
                        atmo_grid[i1].teff,
                    )
        return icor

    def interpolate_corner_models(
        self, teff, logg, monh, icor, atmo_grid, interp="RHOX", itop=1
    ):
        """
        Interpolate over the 8 corner models in the cube around the stellar parameters

        The code below interpolates between 8 corner models to produce
        the output atmosphere. In the first step, pairs of models at each
        of the 4 combinations of log(g) and Teff are interpolated to the
        desired value of [M/H]. These 4 new models are then interpolated
        to the desired value of log(g), yielding 2 models at the requested
        [M/H] and log(g). Finally, this pair of models is interpolated
        to the desired Teff, producing a single output model.

        Interpolation is done on the logarithm of all quantities to improve
        linearity of the fitted data. Kurucz models sometimes have very small
        fractional steps in mass column at the top of the atmosphere. These
        cause wild oscillations in splines fitted to facilitate interpolation
        onto a common depth scale. To circumvent this problem, we ignore the
        top point in the atmosphere by setting itop=1.

        Parameters
        ----------
        teff : float
            effective temperature
        logg : float
            surface gravity
        monh : float
            metallicity
        icor : array_like
            indices of the corner models
        atmo_grid : AtmosphereGrid
            atmosphere grid with input models
        interp : str, optional
            interpolation axis, by default "RHOX"
        itop : int, optional
            index of the topmost point in the RHOX scale, by default 1

        Returns
        -------
        atmo3: Atmosphere
            Interpolated atmosphere
        """

        # We do this for every pair of atmosphere models
        def interpolate(m0, m1, p, param, **kwargs):
            p0 = getattr(m0, param)
            p1 = getattr(m1, param)
            pfrac = (p - p0) / (p1 - p0) if p0 != p1 else 0
            return self.interp_atmo_pair(m0, m1, pfrac, interpvar=interp, **kwargs)

        # Interpolate 8 corner models to create 4 models at the desired [M/H].
        atmo = [[None, None], [None, None]]
        for (i, j) in np.ndindex(2, 2):
            m0 = atmo_grid[icor[0, j, i]]
            m1 = atmo_grid[icor[1, j, i]]
            atmo[i][j] = interpolate(m0, m1, monh, "monh", itop=itop)

        # Interpolate 4 models at the desired [M/H] to create 2 models at desired
        # [M/H] and log(g).
        atmo2 = [None, None]
        for k in range(2):
            atmo2[k] = interpolate(atmo[k][0], atmo[k][1], logg, "logg")

        # Interpolate the 2 models at desired [M/H] and log(g) to create final
        # model at desired [M/H], log(g), and Teff
        atmo3 = interpolate(atmo2[0], atmo2[1], teff, "teff")

        return atmo3

    def spherical_model_correction(self, atmo_grid, icor, logg):
        """
        Correct the radius of the model grids, for the stellar parameters if necessary

        If all interpolated models were spherical, the interpolated model should
        also be reported as spherical. This enables spherical-symmetric radiative
        transfer in the spectral synthesis.

        Formulae for mass and radius at the corners of interpolation cube:
            log(M/Msol) = log g - log g_sol - 2*log(R_sol / R)
            2 * log(R / R_sol) = log g_sol - log g + log(M / M_sol)

        Parameters
        ----------
        atmo_grid : AtmosphereGrid
            complete atmosphere grid
        icor : array_like
            indices for the corner models that were interpolated from the
            atmosphere grid

        Returns
        -------
        geom : {"PP", "SPH"}
            The geometry of the interpolated models
        radius : None, array
            The corrected radius of the models, if geom == "SPH". Otherwise None.
        """

        gtags = atmo_grid.dtype.names
        selection = atmo_grid[icor]
        if "radius" in gtags and "height" in gtags and np.min(selection.radius) > 1:
            mass_cor = (
                selection.logg - logg_sun - 2 * np.log10(R_sun / selection.radius)
            )
            mass = np.mean(mass_cor)
            radius = R_sun * 10 ** ((logg_sun - logg + mass) * 0.5)
            geom = "SPH"
        else:
            geom = "PP"
            radius = None

        return geom, radius

    def interp_atmo_constrained(
        self, x, y, err, par, x2, y2, constraints=None, **kwargs
    ):
        """
        Apply a constraint on each parameter, to have it approach zero

        Parameters
        -------
        x : array[n]
            x data
        y : array[n]
            y data
        err : array[n]
            errors on y data
        par : list[4]
            initial guess for fit parameters - see interp_atmo_func
        x2 : array
            independent variable for tabulated input function
        y2 : array
            dependent variable for tabulated input function
        ** kwargs : dict
            passes keyword arguments to interp_atmo_func

            ndep : int
                number of depth points in supplied quantities
            constraints : array[nconstraint], optional
                error vector for constrained parameters.
                Use errors of 0 for unconstrained parameters.

        Returns
        --------
        ret : list of floats
            best fit parameters
        yfit : array of size (n,)
            best fit to data
        """

        # Evaluate with fixed paramters 3, 4
        # ret = mpfitfun("interp_atmo_func", x, y, err, par, extra=extra, yfit=yfit)
        # TODO: what does the constrain do?
        func = lambda x, *p: self.interp_atmo_func(x, [*p, *par[2:]], x2, y2, **kwargs)
        popt, _ = curve_fit(func, x, y, sigma=err, p0=par[:2])

        ret = [*popt, *par[2:]]
        yfit = func(x, *popt)
        return ret, yfit

    def interp_atmo_func(self, x1, par, x2, y2, ndep=None, y1=None):
        """
        Apply a horizontal shift to x2.
        Apply a vertical shift to y2.
        Interpolate onto x1 the shifted y2 as a function of shifted x2.

        Parameters
        ---------
        x1 : array[ndep1]
            independent variable for output function
        par : array[3]
            shift parameters
            par[0] - horizontal shift for x2
            par[1] - vertical shift for y2
            par[2] - vertical scale factor for y2
        x2 : array[ndep2]
            independent variable for tabulated input function
        y2 : array[ndep2]
            dependent variable for tabulated input function
        ndep : int, optional
            number of depth points in atmospheric structure (default is use all)
        y1 : array[ndep2], optional
            data values being fitted

        Note
        -------
        Only pass y1 if you want to restrict the y-values of extrapolated
        data points. This is useful when solving for the shifts, but should
        not be used when calculating shifted functions for actual use, since
        this restriction can cause discontinuities.
        """
        # Constrained fits may append non-atmospheric quantities to the end of
        # input vector.
        # Extract the output depth scale:
        if ndep is None:
            ndep = len(x1)

        # Shift input x-values.
        # Interpolate input y-values onto output x-values.
        # Shift output y-values.
        x2sh = x2 + par[0]
        y2sh = y2 + par[1]
        y = np.zeros_like(x1)
        y[:ndep] = interp1d(x2sh, y2sh, kind="linear", fill_value="extrapolate")(
            x1[:ndep]
        )  # Note, this implicitly extrapolates

        # Scale output y-values about output y-center.
        ymin = np.min(y[:ndep])
        ymax = np.max(y[:ndep])
        ycen = 0.5 * (ymin + ymax)
        y[:ndep] = ycen + (1.0 + par[2]) * (y[:ndep] - ycen)

        # If extra.y1 was passed, then clip minimum and maximum of output y1.
        if y1 is not None:
            y[:ndep] = np.clip(y[:ndep], np.min(y1), np.max(y1))

        # Set all leftover values to zero
        y[ndep:] = 0
        return y
