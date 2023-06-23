namespace SolPosSharp;


using System.Text;

/// <summary>
/// This class was derived from the NREL Solar Position and Intensity C code:
///
/// https://www.nrel.gov/grid/solar-resource/solpos.html
///
/// Class, function, and variable names were changed to make the code easier to understand, and to
/// follow C# naming convention.
/// </summary>
public static class Solar
{
    /// <summary>
    /// Cumulative number of days prior to beginning of month.
    /// </summary>
    private static readonly int[,] MONTH_DAYS =
    {
        { 0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
        { 0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 },
    };

    /// <summary>
    /// Converts from radians to degrees.
    /// </summary>
    private static readonly double DEG_PER_RAD = 57.295779513;

    /// <summary>
    /// Converts from degrees to radians.
    /// </summary>
    private static readonly double RAD_PER_DEG = 0.0174532925;

    public static SolarPosition CalcSolarPosition(DateTimeOffset time, double latitude, double longitude)
    {
        IoData ioData = new()
        {
            Year = time.Year,
            Month = time.Month,
            DayOfMonth = time.Day,
            Hour = time.Hour,
            Minute = time.Minute,
            Second = time.Second,
            Timezone = time.Offset.TotalHours,
            Latitude = latitude,
            Longitude = longitude,
        };
        // Change the input from day of year to month/day
        ioData.Function &= ~Function.S_DOY;
        ErrorCode result = CalcSolarPosition(ioData);
        if (result == 0)
        {
            return new SolarPosition(ioData.SolarAzimuth, ioData.SolarElevation);
        }
        else
        {
            throw new Exception($"Could not calculate solar position (error code = {result}");
        }
    }

    /// <summary>
    /// Calculates the apparent solar position and the intensity of the sun (theoretical maximum
    /// solar energy) from the time and place on Earth. Inputs are passed in ioData and outputs are
    /// set in ioData.
    ///
    /// Requires (from the struct posdata parameter):
    ///     Date and time:
    ///         year
    ///         daynum   (requirement depends on the S_DOY switch)
    ///         month    (requirement depends on the S_DOY switch)
    ///         day      (requirement depends on the S_DOY switch)
    ///         hour
    ///         minute
    ///         second
    ///         interval  DEFAULT 0
    ///     Location:
    ///         latitude
    ///         longitude
    ///     Location/time adjuster:
    ///         timezone
    ///     Atmospheric pressure and temperature:
    ///         press     DEFAULT 1013.0 mb
    ///         temp      DEFAULT 10.0 degrees C
    ///     Tilt of flat surface that receives solar energy:
    ///         aspect    DEFAULT 180 (South)
    ///         tilt      DEFAULT 0 (Horizontal)
    ///     Function Switch (codes defined in solpos.h)
    ///         function  DEFAULT S_ALL
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <returns>ErrorCodeEnum, which is a bitmask of errors. Zero means no errors.</returns>
    public static ErrorCode CalcSolarPosition(IoData ioData)
    {
        ErrorCode retval;

        TrigData trigData = new TrigData();

        if ((retval = Validate(ioData)) != 0) // validate the inputs
            return retval;

        if (ioData.Function.HasFlag(Function.L_DOY))
            Doy2Dom(ioData); // convert input doy to month-day
        else
            Dom2Doy(ioData); // convert input month-day to doy

        if (ioData.Function.HasFlag(Function.L_GEOM))
            CalcGeometry(ioData); // do basic geometry calculations

        if (ioData.Function.HasFlag(Function.L_ZENETR)) // etr at non-refracted zenith angle
            CalcSolarZenith(ioData, trigData);

        if (ioData.Function.HasFlag(Function.L_SSHA)) // Sunset hour calculation
            CalcSunsetHourAng(ioData, trigData);

        if (ioData.Function.HasFlag(Function.L_SBCF)) // Shadowband correction factor
            CalcShadowbandCorrFact(ioData, trigData);

        if (ioData.Function.HasFlag(Function.L_TST)) // true solar time
            CalcTrueSolarTime(ioData);

        if (ioData.Function.HasFlag(Function.L_SRSS)) // sunrise/sunset calculations
            CalcSunriseSunset(ioData);

        if (ioData.Function.HasFlag(Function.L_SOLAZM)) // solar azimuth calculations
            CalcSolarAzimuth(ioData, trigData);

        if (ioData.Function.HasFlag(Function.L_REFRAC)) // atmospheric refraction calculations
            CalcRefractionCorrection(ioData);

        if (ioData.Function.HasFlag(Function.L_AMASS)) // airmass calculations
            CalcAirMass(ioData);

        if (ioData.Function.HasFlag(Function.L_PRIME)) // kt-prime/unprime calculations
            CalcPrimeUnprime(ioData);

        if (ioData.Function.HasFlag(Function.L_ETR)) // ETR and ETRN (refracted)
            CalcExtTerSolarIrradiance(ioData);

        if (ioData.Function.HasFlag(Function.L_TILT)) // tilt calculations
            CalcExtTerSolIrrTilt(ioData);

        return 0;
    }

    /// <summary>
    /// Validates input parameters.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <returns>ErrorCodeEnum, which is a bitmask of errors. Zero means no errors.</returns>
    private static ErrorCode Validate(IoData ioData)
    {
        ErrorCode retval = 0; // start with no errors

        if (ioData.Function.HasFlag(Function.L_GEOM))
        {
            if (ioData.Year is < 1950 or > 2050) // limits of algorithm
            {
                retval |= ErrorCode.S_YEAR_ERROR;
            }

            if (!ioData.Function.HasFlag(Function.S_DOY) && (ioData.Month < 1 || ioData.Month > 12))
            {
                retval |= ErrorCode.S_MONTH_ERROR;
            }

            if (!ioData.Function.HasFlag(Function.S_DOY) &&
                (ioData.DayOfMonth < 1 || ioData.DayOfMonth > 31))
            {
                retval |= ErrorCode.S_DAY_ERROR;
            }

            if (ioData.Function.HasFlag(Function.S_DOY) &&
                (ioData.DayOfYear < 1 || ioData.DayOfYear > 366))
            {
                retval |= ErrorCode.S_DOY_ERROR;
            }

            if (ioData.Hour is < 0 or > 24)
            {
                retval |= ErrorCode.S_HOUR_ERROR;
            }

            if (ioData.Minute is < 0 or > 59)
            {
                retval |= ErrorCode.S_MINUTE_ERROR;
            }

            if (ioData.Second is < 0 or > 59)
            {
                retval |= ErrorCode.S_SECOND_ERROR;
            }

            if (ioData.Hour == 24 && ioData.Minute > 0) // no more than 24 hrs
                retval |= ErrorCode.S_HOUR_ERROR | ErrorCode.S_MINUTE_ERROR;

            if (ioData.Hour == 24 && ioData.Second > 0) // no more than 24 hrs
                retval |= ErrorCode.S_HOUR_ERROR | ErrorCode.S_SECOND_ERROR;

            if (Math.Abs(ioData.Timezone) > 12.0)
            {
                retval |= ErrorCode.S_TZONE_ERROR;
            }

            if (ioData.Interval is < 0 or > 28800)
            {
                retval |= ErrorCode.S_INTRVL_ERROR;
            }

            if (Math.Abs(ioData.Longitude) > 180.0)
            {
                retval |= ErrorCode.S_LON_ERROR;
            }

            if (Math.Abs(ioData.Latitude) > 90.0)
            {
                retval |= ErrorCode.S_LAT_ERROR;
            }
        }

        if (ioData.Function.HasFlag(Function.L_REFRAC) && Math.Abs(ioData.Temp) > 100.0)
        {
            retval |= ErrorCode.S_TEMP_ERROR;
        }

        if ((ioData.Function.HasFlag(Function.L_REFRAC) && ioData.Press < 0.0) ||
            ioData.Press > 2000.0)
        {
            retval |= ErrorCode.S_PRESS_ERROR;
        }

        if (ioData.Function.HasFlag(Function.L_TILT) && Math.Abs(ioData.PanelTilt) > 180.0)
        {
            retval |= ErrorCode.S_TILT_ERROR;
        }

        if (ioData.Function.HasFlag(Function.L_TILT) && Math.Abs(ioData.PanelAzimuth) > 360.0)
        {
            retval |= ErrorCode.S_ASPECT_ERROR;
        }

        if ((ioData.Function.HasFlag(Function.L_SBCF) && ioData.ShadowBandWidth < 1.0) ||
            ioData.ShadowBandWidth > 100.0)
        {
            retval |= ErrorCode.S_SBWID_ERROR;
        }

        if ((ioData.Function.HasFlag(Function.L_SBCF) && ioData.ShadowBandRadius < 1.0) ||
            ioData.ShadowBandRadius > 100.0)
        {
            retval |= ErrorCode.S_SBRAD_ERROR;
        }

        if (ioData.Function.HasFlag(Function.L_SBCF) && Math.Abs(ioData.ShadowBandSkyFact) > 1.0)
        {
            retval |= ErrorCode.S_SBSKY_ERROR;
        }

        return retval;
    }

    /// <summary>
    /// Calculate the month and day of month from the day of year.
    /// 
    /// Requires (from struct posdata parameter):
    ///     Year and day number:
    ///         year
    ///         daynum
    ///
    /// Returns (via the struct posdata parameter):
    ///         year
    ///         month
    ///         day
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void Doy2Dom(IoData ioData)
    {
        int imon; // Month (month_days) array counter
        int leap; // leap year switch

        /* Set the leap year switch */
        if (((ioData.Year % 4) == 0) && (((ioData.Year % 100) != 0) || ((ioData.Year % 400) == 0)))
            leap = 1;
        else
            leap = 0;

        /* Find the month */
        imon = 12;
        while (ioData.DayOfYear <= MONTH_DAYS[leap, imon])
        {
            --imon;
        }

        /* Set the month and day of month */
        ioData.Month = imon;
        ioData.DayOfMonth = ioData.DayOfYear - MONTH_DAYS[leap, imon];
    }

    /// <summary>
    /// Calculates day of year from day of month.
    /// 
    /// Requires (from struct posdata parameter):
    ///         year
    ///         month
    ///         day
    ///
    /// Returns (via the struct posdata parameter):
    ///         year
    ///         daynum
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void Dom2Doy(IoData ioData)
    {
        ioData.DayOfYear = ioData.DayOfMonth + MONTH_DAYS[0, ioData.Month];

        /* (adjust for leap year) */
        if (ioData.Year % 4 == 0 && (ioData.Year % 100 != 0 || ioData.Year % 400 == 0) &&
            ioData.Month > 2)
            ioData.DayOfYear += 1;
    }

    /// <summary>
    /// Calculates underlying geometry for a given time and location.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcGeometry(IoData ioData)
    {
        double denom; // denominator (bottom) of the fraction
        double cosD2; // cosine of d2
        double cosDayAng; // cosine of the day angle or delination
        double cosDayAng2; // pdat->dayang times two
        double deltaYear; // difference between current year and 1949
        double sinD2; // sine of d2
        double sinDayAng; // sine of the day angle
        double numer; // numerator (top) of the fraction
        int leapYrCntr; // leap year counter

        /* Day angle */
        /*  Iqbal, M.  1983.  An Introduction to Solar Radiation.
              Academic Press, NY., page 3 */
        ioData.DayAngle = 360.0 * (ioData.DayOfYear - 1) / 365.0;

        /* Earth radius vector * solar constant = solar energy */
        /*  Spencer, J. W.  1971.  Fourier series representation of the
            position of the sun.  Search 2 (5), page 172 */
        sinDayAng = Math.Sin(RAD_PER_DEG * ioData.DayAngle);
        cosDayAng = Math.Cos(RAD_PER_DEG * ioData.DayAngle);
        cosDayAng2 = 2.0 * ioData.DayAngle;
        cosD2 = Math.Cos(RAD_PER_DEG * cosDayAng2);
        sinD2 = Math.Sin(RAD_PER_DEG * cosDayAng2);

        ioData.EarthRadiusVector = 1.000110 + (0.034221 * cosDayAng) + (0.001280 * sinDayAng);
        ioData.EarthRadiusVector += (0.000719 * cosD2) + (0.000077 * sinD2);

        /* Universal Coordinated (Greenwich standard) time */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.TimeUtc = (ioData.Hour * 3600.0) + (ioData.Minute * 60.0) + ioData.Second -
                         (ioData.Interval / 2.0);
        ioData.TimeUtc = (ioData.TimeUtc / 3600.0) - ioData.Timezone;

        /* Julian Day minus 2,400,000 days (to eliminate roundoff errors) */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */

        /* No adjustment for century non-leap years since this function is
           bounded by 1950 - 2050 */
        deltaYear = ioData.Year - 1949;
        leapYrCntr = (int)(deltaYear / 4.0);
        ioData.JulianDay = 32916.5 + (deltaYear * 365.0) + leapYrCntr + ioData.DayOfYear +
                           (ioData.TimeUtc / 24.0);

        /* Time used in the calculation of ecliptic coordinates */
        /* Noon 1 JAN 2000 = 2,400,000 + 51,545 days Julian Date */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.EclipticTime = ioData.JulianDay - 51545.0;

        /* Mean longitude */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.MeanLon = 280.460 + (0.9856474 * ioData.EclipticTime);

        /* (dump the multiples of 360, so the answer is between 0 and 360) */
        ioData.MeanLon -= 360.0 * (int)(ioData.MeanLon / 360.0);
        if (ioData.MeanLon < 0.0)
            ioData.MeanLon += 360.0;

        /* Mean anomaly */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.MeanAnom = 357.528 + (0.9856003 * ioData.EclipticTime);

        /* (dump the multiples of 360, so the answer is between 0 and 360) */
        ioData.MeanAnom -= 360.0 * (int)(ioData.MeanAnom / 360.0);
        if (ioData.MeanAnom < 0.0)
            ioData.MeanAnom += 360.0;

        /* Ecliptic longitude */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.EclipticLon = ioData.MeanLon + (1.915 * Math.Sin(ioData.MeanAnom * RAD_PER_DEG)) +
                             (0.020 * Math.Sin(2.0 * ioData.MeanAnom * RAD_PER_DEG));

        /* (dump the multiples of 360, so the answer is between 0 and 360) */
        ioData.EclipticLon -= 360.0 * (int)(ioData.EclipticLon / 360.0);
        if (ioData.EclipticLon < 0.0)
            ioData.EclipticLon += 360.0;

        /* Obliquity of the ecliptic */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */

        /* 02 Feb 2001 SMW corrected sign in the following line */
        /*  pdat->ecobli = 23.439 + 4.0e-07 * pdat->ectime;     */
        ioData.EclipticObliquity = 23.439 - (4.0e-07 * ioData.EclipticTime);

        /* Declination */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.Declination = DEG_PER_RAD *
                             Math.Asin(Math.Sin(ioData.EclipticObliquity * RAD_PER_DEG) *
                                       Math.Sin(ioData.EclipticLon * RAD_PER_DEG));

        /* Right ascension */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        numer = Math.Cos(RAD_PER_DEG * ioData.EclipticObliquity) *
                Math.Sin(RAD_PER_DEG * ioData.EclipticLon);
        denom = Math.Cos(RAD_PER_DEG * ioData.EclipticLon);

        ioData.RightAscension = DEG_PER_RAD * Math.Atan2(numer, denom);

        /* (make it a positive angle) */
        if (ioData.RightAscension < 0.0)
            ioData.RightAscension += 360.0;

        /* Greenwich mean sidereal time */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.Gmst = 6.697375 + (0.0657098242 * ioData.EclipticTime) + ioData.TimeUtc;

        /* (dump the multiples of 24, so the answer is between 0 and 24) */
        ioData.Gmst -= 24.0 * (int)(ioData.Gmst / 24.0);
        if (ioData.Gmst < 0.0)
            ioData.Gmst += 24.0;

        /* Local mean sidereal time */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.Lmst = (ioData.Gmst * 15.0) + ioData.Longitude;

        /* (dump the multiples of 360, so the answer is between 0 and 360) */
        ioData.Lmst -= 360.0 * (int)(ioData.Lmst / 360.0);
        if (ioData.Lmst < 0.0)
            ioData.Lmst += 360.0;

        /* Hour angle */
        /*  Michalsky, J.  1988.  The Astronomical Almanac's algorithm for
            approximate solar position (1950-2050).  Solar Energy 40 (3),
            pp. 227-235. */
        ioData.HourAng = ioData.Lmst - ioData.RightAscension;

        /* (force it between -180 and 180 degrees) */
        if (ioData.HourAng < -180.0)
            ioData.HourAng += 360.0;
        else if (ioData.HourAng > 180.0)
            ioData.HourAng -= 360.0;
    }

    /// <summary>
    /// Calculate (uncorrected) solar zenith angle.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <param name="trigData">Structure used for trigonometric data.</param>
    private static void CalcSolarZenith(IoData ioData, TrigData trigData)
    {
        double cz; // cosine of the solar zenith angle

        CalcLocalTrig(ioData, trigData);
        cz = (trigData.SinDeclination * trigData.SinLat) +
             (trigData.CosDeclination * trigData.CosLat * trigData.CosHourAng);

        /* (watch out for the roundoff errors) */
        if (Math.Abs(cz) > 1.0)
        {
            if (cz >= 0.0)
                cz = 1.0;
            else
                cz = -1.0;
        }

        ioData.SolarZenith = Math.Acos(cz) * DEG_PER_RAD;

        /* (limit the degrees below the horizon to 9 [+90 -> 99]) */
        if (ioData.SolarZenith > 99.0)
            ioData.SolarZenith = 99.0;

        ioData.SolarElevation = 90.0 - ioData.SolarZenith;
    }

    /// <summary>
    /// Calculate trigonometric values used by several functions.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <param name="trigData">Structure used for trigonometric data.</param>
    private static void CalcLocalTrig(IoData ioData, TrigData trigData)
    {
        if (trigData.SinDeclination < -900.0) // sd was initialized -999 as flag
        {
            trigData.SinDeclination = 1.0; // reflag as having completed calculations
            if ((ioData.Function | Function.CD_MASK) != 0)
                trigData.CosDeclination = Math.Cos(RAD_PER_DEG * ioData.Declination);

            if ((ioData.Function | Function.CH_MASK) != 0)
                trigData.CosHourAng = Math.Cos(RAD_PER_DEG * ioData.HourAng);

            if ((ioData.Function | Function.CL_MASK) != 0)
                trigData.CosLat = Math.Cos(RAD_PER_DEG * ioData.Latitude);

            if ((ioData.Function | Function.SD_MASK) != 0)
                trigData.SinDeclination = Math.Sin(RAD_PER_DEG * ioData.Declination);

            if ((ioData.Function | Function.SL_MASK) != 0)
                trigData.SinLat = Math.Sin(RAD_PER_DEG * ioData.Latitude);
        }
    }

    /// <summary>
    /// Calculate sunset hour angle in degrees.
    /// 
    /// Iqbal, M.  1983.  An Introduction to Solar Radiation.
    ///     Academic Press, NY., page 16
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <param name="trigData">Structure used for trigonometric data.</param>
    private static void CalcSunsetHourAng(IoData ioData, TrigData trigData)
    {
        double cssha; // cosine of the sunset hour angle
        double cdcl; // ( cd * cl )

        CalcLocalTrig(ioData, trigData);
        cdcl = trigData.CosDeclination * trigData.CosLat;

        if (Math.Abs(cdcl) >= 0.001)
        {
            cssha = -trigData.SinLat * trigData.SinDeclination / cdcl;

            /* This keeps the cosine from blowing on roundoff */
            if (cssha < -1.0)
                ioData.SunsetHourAng = 180.0;
            else if (cssha > 1.0)
                ioData.SunsetHourAng = 0.0;
            else
                ioData.SunsetHourAng = DEG_PER_RAD * Math.Acos(cssha);
        }
        else if (((ioData.Declination >= 0.0) && (ioData.Latitude > 0.0)) ||
                 ((ioData.Declination < 0.0) && (ioData.Latitude < 0.0)))
            ioData.SunsetHourAng = 180.0;
        else
            ioData.SunsetHourAng = 0.0;
    }

    /// <summary>
    /// Calculate shadowband correction factor.
    /// 
    /// Drummond, A. J.  1956.  A contribution to absolute pyrheliometry.
    ///     Q. J. R. Meteorol. Soc. 82, pp. 481-493
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <param name="trigData">Structure used for trigonometric data.</param>
    private static void CalcShadowbandCorrFact(IoData ioData, TrigData trigData)
    {
        double p; // used to compute sbcf
        double t1;
        double t2;

        CalcLocalTrig(ioData, trigData);
        p = 0.6366198 * ioData.ShadowBandWidth / ioData.ShadowBandRadius *
            Math.Pow(trigData.CosDeclination, 3);
        t1 = trigData.SinLat * trigData.SinDeclination * ioData.SunsetHourAng * RAD_PER_DEG;
        t2 = trigData.CosLat * trigData.CosDeclination *
             Math.Sin(ioData.SunsetHourAng * RAD_PER_DEG);
        ioData.ShadowBandCorrFact = ioData.ShadowBandSkyFact + (1.0 / (1.0 - (p * (t1 + t2))));
    }

    /// <summary>
    /// Calculate true solar time (TST) in minutes from midnight.
    /// 
    /// Iqbal, M.  1983.  An Introduction to Solar Radiation.
    ///     Academic Press, NY., page 13
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcTrueSolarTime(IoData ioData)
    {
        ioData.TrueSolarTime = (180.0 + ioData.HourAng) * 4.0;
        ioData.TrueSolarTimeLoc =
            ioData.TrueSolarTime -
            (ioData.Hour * 60.0) -
            ioData.Minute -
            (ioData.Second / 60.0) +
            (ioData.Interval / 120.0); // add back half of the interval

        /* bound tstfix to this day */
        while (ioData.TrueSolarTimeLoc > 720.0)
            ioData.TrueSolarTimeLoc -= 1440.0;

        while (ioData.TrueSolarTimeLoc < -720.0)
            ioData.TrueSolarTimeLoc += 1440.0;

        ioData.EqOfTime = ioData.TrueSolarTimeLoc + (60.0 * ioData.Timezone) -
                          (4.0 * ioData.Longitude);
    }

    /// <summary>
    /// Calculate sunrise and sunset times in minutes from midnight.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcSunriseSunset(IoData ioData)
    {
        if (ioData.SunsetHourAng <= 1.0)
        {
            ioData.SunriseTime = 2999.0;
            ioData.SunsetTime = -2999.0;
        }
        else if (ioData.SunsetHourAng >= 179.0)
        {
            ioData.SunriseTime = -2999.0;
            ioData.SunsetTime = 2999.0;
        }
        else
        {
            ioData.SunriseTime = 720.0 - (4.0 * ioData.SunsetHourAng) - ioData.TrueSolarTimeLoc;
            ioData.SunsetTime = 720.0 + (4.0 * ioData.SunsetHourAng) - ioData.TrueSolarTimeLoc;
        }
    }

    /// <summary>
    /// Calculate solar azimuth angle in degrees.
    /// 
    /// Iqbal, M.  1983.  An Introduction to Solar Radiation.
    ///     Academic Press, NY., page 13
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <param name="trigData">Structure used for trigonometric data.</param>
    private static void CalcSolarAzimuth(IoData ioData, TrigData trigData)
    {
        double ca; // cosine of the solar azimuth angle
        double ce; // cosine of the solar elevation
        double cecl; // ( ce * cl )
        double se; // sine of the solar elevation

        CalcLocalTrig(ioData, trigData);
        ce = Math.Cos(RAD_PER_DEG * ioData.SolarElevation);
        se = Math.Sin(RAD_PER_DEG * ioData.SolarElevation);

        ioData.SolarAzimuth = 180.0;
        cecl = ce * trigData.CosLat;
        if (Math.Abs(cecl) >= 0.001)
        {
            ca = ((se * trigData.SinLat) - trigData.SinDeclination) / cecl;
            if (ca > 1.0)
            {
                ca = 1.0;
            }
            else if (ca < -1.0)
            {
                ca = -1.0;
            }

            ioData.SolarAzimuth = 180.0 - (Math.Acos(ca) * DEG_PER_RAD);
            if (ioData.HourAng > 0)
            {
                ioData.SolarAzimuth = 360.0 - ioData.SolarAzimuth;
            }
        }
    }

    /// <summary>
    /// Calculate refraction correction in degrees.
    /// 
    /// Zimmerman, John C.  1981.  Sun-pointing programs and their accuracy.
    ///     SAND81-0761, Experimental Systems Operation Division 4721,
    ///     Sandia National Laboratories, Albuquerque, NM.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcRefractionCorrection(IoData ioData)
    {
        double prestemp; // temporary pressure/temperature correction
        double refcor; // temporary refraction correction
        double tanelev; // tangent of the solar elevation angle

        /* If the sun is near zenith, the algorithm bombs; refraction near 0 */
        if (ioData.SolarElevation > 85.0)
            refcor = 0.0;
        /* Otherwise, we have refraction */
        else
        {
            tanelev = Math.Tan(RAD_PER_DEG * ioData.SolarElevation);
            if (ioData.SolarElevation >= 5.0)
            {
                refcor = (58.1 / tanelev) - (0.07 / Math.Pow(tanelev, 3)) +
                         (0.000086 / Math.Pow(tanelev, 5));
            }
            else if (ioData.SolarElevation >= -0.575)
            {
                refcor = 1735.0 + (ioData.SolarElevation * (-518.2 + (ioData.SolarElevation *
                    (103.4 + (ioData.SolarElevation *
                              (-12.79 + (ioData.SolarElevation * 0.711)))))));
            }
            else
            {
                refcor = -20.774 / tanelev;
            }

            prestemp = ioData.Press * 283.0 / (1013.0 * (273.0 + ioData.Temp));
            refcor *= prestemp / 3600.0;
        }

        /* Refracted solar elevation angle */
        ioData.SolarElevationCorr = ioData.SolarElevation + refcor;

        /* (limit the degrees below the horizon to 9) */
        if (ioData.SolarElevationCorr < -9.0)
        {
            ioData.SolarElevationCorr = -9.0;
        }

        /* Refracted solar zenith angle */
        ioData.SolarZenithCorr = 90.0 - ioData.SolarElevationCorr;
        ioData.CosSolarZenith = Math.Cos(RAD_PER_DEG * ioData.SolarZenithCorr);
    }

    /// <summary>
    /// Calculate air mass.
    /// 
    /// Kasten, F. and Young, A.  1989.  Revised optical air mass tables and
    ///     approximation formula.  Applied Optics 28 (22), pp. 4735-4738.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcAirMass(IoData ioData)
    {
        if (ioData.SolarZenithCorr > 93.0)
        {
            ioData.AirMass = -1.0;
            ioData.AirMassPress = -1.0;
        }
        else
        {
            ioData.AirMass = 1.0 / (Math.Cos(RAD_PER_DEG * ioData.SolarZenithCorr) +
                                    (0.50572 * Math.Pow(96.07995 - ioData.SolarZenithCorr,
                                        -1.6364)));
            ioData.AirMassPress = ioData.AirMass * ioData.Press / 1013.0;
        }
    }

    /// <summary>
    /// Calculate prime factor used to convert Kt to normalized Kt', and unprime factor, used to
    /// convert Kt' to Kt.
    /// 
    /// Perez, R., P. Ineichen, Seals, R., & Zelenka, A.  1990.  Making full use of the clearness
    ///     index for parameterizing hourly insolation conditions. Solar Energy 45 (2), pp. 111-114
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcPrimeUnprime(IoData ioData)
    {
        ioData.Unprime = (1.031 * Math.Exp(-1.4 / (0.9 + (9.4 / ioData.AirMass)))) + 0.1;
        ioData.Prime = 1.0 / ioData.Unprime;
    }

    /// <summary>
    /// Calculate extraterrestrial (top-of-atmosphere) solar irradiance.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcExtTerSolarIrradiance(IoData ioData)
    {
        if (ioData.CosSolarZenith > 0.0)
        {
            ioData.ExtTerDNI = ioData.SolarConst * ioData.EarthRadiusVector;
            ioData.ExtTerGHI = ioData.ExtTerDNI * ioData.CosSolarZenith;
        }
        else
        {
            ioData.ExtTerDNI = 0.0;
            ioData.ExtTerGHI = 0.0;
        }
    }

    /// <summary>
    /// Calculate extraterrestrial (top-of-atmosphere) solar irradiance on tilted surface.
    /// </summary>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    private static void CalcExtTerSolIrrTilt(IoData ioData)
    {
        double ca; // cosine of the solar azimuth angle
        double cp; // cosine of the panel aspect
        double ct; // cosine of the panel tilt
        double sa; // sine of the solar azimuth angle
        double sp; // sine of the panel aspect
        double st; // sine of the panel tilt
        double sz; // sine of the refraction corrected solar zenith angle

        /* Cosine of the angle between the sun and a tipped flat surface,
           useful for calculating solar energy on tilted surfaces */
        ca = Math.Cos(RAD_PER_DEG * ioData.SolarAzimuth);
        cp = Math.Cos(RAD_PER_DEG * ioData.PanelAzimuth);
        ct = Math.Cos(RAD_PER_DEG * ioData.PanelTilt);
        sa = Math.Sin(RAD_PER_DEG * ioData.SolarAzimuth);
        sp = Math.Sin(RAD_PER_DEG * ioData.PanelAzimuth);
        st = Math.Sin(RAD_PER_DEG * ioData.PanelTilt);
        sz = Math.Sin(RAD_PER_DEG * ioData.SolarZenithCorr);
        ioData.CosSolarIncidence =
            (ioData.CosSolarZenith * ct) + (sz * st * ((ca * cp) + (sa * sp)));

        if (ioData.CosSolarIncidence > 0.0)
        {
            ioData.ExtTerTilt = ioData.ExtTerDNI * ioData.CosSolarIncidence;
        }
        else
        {
            ioData.ExtTerTilt = 0.0;
        }
    }

    /// <summary>
    /// Decode error code into friendly messages.
    /// </summary>
    /// <param name="code">Error code to decode</param>
    /// <param name="ioData">Structure used for inputs, outputs, and transitional values.</param>
    /// <returns>String containing error messages.</returns>
    public static string DecodeErrorCode(ErrorCode code, IoData ioData)
    {
        StringBuilder sb = new StringBuilder();

        if (code.HasFlag(ErrorCode.S_YEAR_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the year: {ioData.Year:D} [1950-2050]\n");

        if (code.HasFlag(ErrorCode.S_MONTH_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the month: {ioData.Month:D}\n");

        if (code.HasFlag(ErrorCode.S_DAY_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the day-of-month: {ioData.DayOfMonth:D}\n");

        if (code.HasFlag(ErrorCode.S_DOY_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the day-of-year: {ioData.DayOfYear:D}\n");

        if (code.HasFlag(ErrorCode.S_HOUR_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the hour: {ioData.Hour:D}\n");

        if (code.HasFlag(ErrorCode.S_MINUTE_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the minute: {ioData.Minute:D}\n");

        if (code.HasFlag(ErrorCode.S_SECOND_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the second: {ioData.Second:D}\n");

        if (code.HasFlag(ErrorCode.S_TZONE_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the time zone: {ioData.Timezone:f}\n");

        if (code.HasFlag(ErrorCode.S_INTRVL_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the interval: {ioData.Interval:D}\n");

        if (code.HasFlag(ErrorCode.S_LAT_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the latitude: {ioData.Latitude:f}\n");

        if (code.HasFlag(ErrorCode.S_LON_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the longitude: {ioData.Longitude:f}\n");

        if (code.HasFlag(ErrorCode.S_TEMP_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the temperature: {ioData.Temp:f}\n");

        if (code.HasFlag(ErrorCode.S_PRESS_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the pressure: {ioData.Press:f}\n");

        if (code.HasFlag(ErrorCode.S_TILT_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the tilt: {ioData.PanelTilt:f}\n");

        if (code.HasFlag(ErrorCode.S_ASPECT_ERROR))
            sb.AppendLine($"S_decode ==> Please fix the aspect: {ioData.PanelAzimuth:f}\n");

        if (code.HasFlag(ErrorCode.S_SBWID_ERROR))
            sb.AppendLine(
                $"S_decode ==> Please fix the shadowband width: {ioData.ShadowBandWidth:f}\n");

        if (code.HasFlag(ErrorCode.S_SBRAD_ERROR))
            sb.AppendLine(
                $"S_decode ==> Please fix the shadowband radius: {ioData.ShadowBandRadius:f}\n");

        if (code.HasFlag(ErrorCode.S_SBSKY_ERROR))
            sb.AppendLine(
                $"S_decode ==> Please fix the shadowband sky factor: {ioData.ShadowBandSkyFact:f}\n");

        return sb.ToString();
    }
}

[Flags]
public enum ErrorCode
{
    /*    Code                  Bit        Parameter            Range
    =============             ======= ===================   =============   */
    S_YEAR_ERROR = 1 << 0, //  0   year                  1950 -  2050
    S_MONTH_ERROR = 1 << 1, //  1   month                    1 -    12
    S_DAY_ERROR = 1 << 2, //  2   day-of-month             1 -    31
    S_DOY_ERROR = 1 << 3, //  3   day-of-year              1 -   366
    S_HOUR_ERROR = 1 << 4, //  4   hour                     0 -    24
    S_MINUTE_ERROR = 1 << 5, //  5   minute                   0 -    59
    S_SECOND_ERROR = 1 << 6, //  6   second                   0 -    59
    S_TZONE_ERROR = 1 << 7, //  7   time zone              -12 -    12
    S_INTRVL_ERROR = 1 << 8, //  8   interval (seconds)       0 - 28800
    S_LAT_ERROR = 1 << 9, //  9   latitude               -90 -    90
    S_LON_ERROR = 1 << 10, // 10   longitude             -180 -   180
    S_TEMP_ERROR = 1 << 11, // 11   temperature (deg. C)  -100 -   100
    S_PRESS_ERROR = 1 << 12, // 12   pressure (millibars)     0 -  2000
    S_TILT_ERROR = 1 << 13, // 13   tilt                   -90 -    90
    S_ASPECT_ERROR = 1 << 14, // 14   aspect                -360 -   360
    S_SBWID_ERROR = 1 << 15, // 15   shadow band width (cm)   1 -   100
    S_SBRAD_ERROR = 1 << 16, // 16   shadow band radius (cm)  1 -   100
    S_SBSKY_ERROR = 1 << 17 // 17   shadow band sky factor  -1 -     1
}

[Flags]
public enum Function
{
    L_DOY = 1 << 1, // L_DOY = 0x0001;
    L_GEOM = 1 << 2, // L_GEOM = 0x0002;
    L_ZENETR = 1 << 3, // L_ZENETR = 0x0004;
    L_SSHA = 1 << 4, // L_SSHA = 0x0008;
    L_SBCF = 1 << 5, // L_SBCF = 0x0010;
    L_TST = 1 << 6, // L_TST = 0x0020;
    L_SRSS = 1 << 7, // L_SRSS = 0x0040;
    L_SOLAZM = 1 << 8, // L_SOLAZM = 0x0080;
    L_REFRAC = 1 << 9, // L_REFRAC = 0x0100; 
    L_AMASS = 1 << 10, // L_AMASS = 0x0200;
    L_PRIME = 1 << 11, // L_PRIME = 0x0400;
    L_TILT = 1 << 12, // L_TILT = 0x0800;
    L_ETR = 1 << 13, // L_ETR = 0x1000;

    L_DEFAULT = L_GEOM | L_ZENETR | L_SSHA | L_SBCF | L_TST | L_SRSS | L_SOLAZM | L_REFRAC |
                L_AMASS | L_PRIME | L_TILT | L_ETR,

    L_ALL = L_DOY | L_GEOM | L_ZENETR | L_SSHA | L_SBCF | L_TST | L_SRSS | L_SOLAZM | L_REFRAC |
            L_AMASS | L_PRIME | L_TILT | L_ETR, // L_ALL = 0xFFFF;

    S_DOY = L_DOY,
    S_GEOM = L_GEOM | S_DOY,
    S_ZENETR = L_ZENETR | S_GEOM,
    S_SSHA = L_SSHA | S_GEOM,
    S_SBCF = L_SBCF | S_SSHA,
    S_TST = L_TST | S_GEOM,
    S_SRSS = L_SRSS | S_SSHA | S_TST,
    S_SOLAZM = L_SOLAZM | S_ZENETR,
    S_REFRAC = L_REFRAC | S_ZENETR,
    S_AMASS = L_AMASS | S_REFRAC,
    S_PRIME = L_PRIME | S_AMASS,
    S_TILT = L_TILT | S_SOLAZM | S_REFRAC,
    S_ETR = L_ETR | S_REFRAC,
    S_ALL = L_ALL,

    SD_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
    SL_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
    CL_MASK = L_ZENETR | L_SSHA | S_SBCF | S_SOLAZM,
    CD_MASK = L_ZENETR | L_SSHA | S_SBCF,
    CH_MASK = L_ZENETR
}

public class TrigData // used to pass calculated values locally
{
    public TrigData()
    {
        /* initialize the trig structure */
        SinDeclination = -999.0; // flag to force calculation of trig data
        CosDeclination = 1.0;
        CosHourAng = 1.0; // set the rest of these to something safe
        CosLat = 1.0;
        SinLat = 1.0;
    }

    public double CosDeclination { get; set; } // cosine of the declination
    public double CosHourAng { get; set; } // cosine of the hour angle
    public double CosLat { get; set; } // cosine of the latitude
    public double SinDeclination { get; set; } // sine of the declination
    public double SinLat { get; set; } // sine of the latitude
}

/// <summary>
/// Structure used for inputs, outputs, and transitional variables. It is passed around by
/// reference.
/// </summary>
public class IoData
{
    public IoData()
    {
        DayOfMonth = -99;
        DayOfYear = -999;
        Hour = -99;
        Minute = -99;
        Month = -99;
        Second = -99;
        Year = -99;
        Interval = 0;
        PanelAzimuth = 180.0;
        Latitude = -99.0;
        Longitude = -999.0;
        Press = 1013.0;
        SolarConst = 1367.0;
        Temp = 15.0;
        PanelTilt = 0.0;
        Timezone = -99.0;
        ShadowBandWidth = 7.6;
        ShadowBandRadius = 31.7;
        ShadowBandSkyFact = 0.04;
        Function = Function.S_ALL;
    }

    /// <summary>
    /// I/O: S_DOY Day of month (May 27 = 27, etc.) solpos will CALCULATE this by default, or will 
    /// optionally require it as input depending on the setting of the S_DOY function switch.
    /// </summary>
    public int DayOfMonth { get; set; }

    /// <summary>
    /// I/O: S_DOY Day number (day of year; Feb 1 = 32 ) solpos REQUIRES this by default, but will
    /// optionally calculate it from month and day depending on the setting of the S_DOY function
    /// switch.
    /// </summary>
    public int DayOfYear { get; set; }

    /// <summary>
    /// I: Bitmap to choose functions for desired output.
    /// </summary>
    public Function Function { get; set; }

    /// <summary>
    /// I: Hour of day, 0 - 23, DEFAULT = 12
    /// </summary>
    public int Hour { get; set; }

    /// <summary>
    /// I: Interval of a measurement period in seconds. Forces solpos to use the time and date from
    /// the interval midpoint. The INPUT time (hour, minute, and second) is assumed to be the END
    /// of the measurement interval.
    /// </summary>
    public int Interval { get; set; }

    /// <summary>
    /// I: Minute of hour, 0 - 59, DEFAULT = 0
    /// </summary>
    public int Minute { get; set; }

    /// <summary>
    /// I/O: S_DOY Month number (Jan = 1, Feb = 2, etc.) solpos will CALCULATE this by default, or
    /// will optionally require it as input depending on the setting of the S_DOY function switch.
    /// </summary>
    public int Month { get; set; }

    /// <summary>
    /// I: Second of minute, 0 - 59, DEFAULT = 0
    /// </summary>
    public int Second { get; set; }

    /// <summary>
    /// I: 4-digit year (2-digit year is NOT allowed
    /// </summary>
    public int Year { get; set; }

    /// <summary>
    /// O: S_AMASS Relative optical airmass
    /// </summary>
    public double AirMass { get; set; }

    /// <summary>
    /// O: S_AMASS Pressure-corrected airmass
    /// </summary>
    public double AirMassPress { get; set; }

    /// <summary>
    /// I: Azimuth of panel surface (direction it faces) N=0, E=90, S=180, W=270, DEFAULT = 180
    /// </summary>
    public double PanelAzimuth { get; set; }

    /// <summary>
    /// O: S_SOLAZM Solar azimuth angle:  N=0, E=90, S=180, W=270
    /// </summary>
    public double SolarAzimuth { get; set; }

    /// <summary>
    /// O: S_TILT Cosine of solar incidence angle on panel
    /// </summary>
    public double CosSolarIncidence { get; set; }

    /// <summary>
    /// O: S_REFRAC Cosine of refraction corrected solar zenith angle
    /// </summary>
    public double CosSolarZenith { get; set; }

    /// <summary>
    /// T: S_GEOM Day angle (daynum*360/year-length) degrees
    /// </summary>
    public double DayAngle { get; set; }

    /// <summary>
    /// T: S_GEOM Declination--zenith angle of solar noon at equator, degrees NORTH
    /// </summary>
    public double Declination { get; set; }

    /// <summary>
    /// T: S_GEOM Ecliptic longitude, degrees
    /// </summary>
    public double EclipticLon { get; set; }

    /// <summary>
    /// T: S_GEOM Obliquity of ecliptic
    /// </summary>
    public double EclipticObliquity { get; set; }

    /// <summary>
    /// T: S_GEOM Time of ecliptic calculations
    /// </summary>
    public double EclipticTime { get; set; }

    /// <summary>
    /// O: S_ZENETR Solar elevation, no atmospheric correction (= ETR)
    /// </summary>
    public double SolarElevation { get; set; }

    /// <summary>
    /// O: S_REFRAC Solar elevation angle, deg. from horizon, refracted
    /// </summary>
    public double SolarElevationCorr { get; set; }

    /// <summary>
    /// T: S_TST Equation of time (TST - LMT), minutes
    /// </summary>
    public double EqOfTime { get; set; }

    /// <summary>
    /// T: S_GEOM Earth radius vector (multiplied to solar constant)
    /// </summary>
    public double EarthRadiusVector { get; set; }

    /// <summary>
    /// O: S_ETR Extraterrestrial (top-of-atmosphere) W/sq m global horizontal solar irradiance
    /// </summary>
    public double ExtTerGHI { get; set; }

    /// <summary>
    /// O: S_ETR Extraterrestrial (top-of-atmosphere) W/sq m direct normal solar  irradiance
    /// </summary>
    public double ExtTerDNI { get; set; }

    /// <summary>
    /// O: S_TILT Extraterrestrial (top-of-atmosphere) W/sq m global irradiance on a tilted surface
    /// </summary>
    public double ExtTerTilt { get; set; }

    /// <summary>
    /// T: S_GEOM Greenwich mean sidereal time, hours
    /// </summary>
    public double Gmst { get; set; }

    /// <summary>
    /// T: S_GEOM Hour angle--hour of sun from solar noon, degrees WEST
    /// </summary>
    public double HourAng { get; set; }

    /// <summary>
    /// T: S_GEOM Julian Day of 1 JAN 2000 minus 2,400,000 days (in order to regain single precision)
    /// </summary>
    public double JulianDay { get; set; }

    /// <summary>
    /// I: Latitude, degrees north (south negative)
    /// </summary>
    public double Latitude { get; set; }

    /// <summary>
    /// I: Longitude, degrees east (west negative)
    /// </summary>
    public double Longitude { get; set; }

    /// <summary>
    /// T: S_GEOM Local mean sidereal time, degrees
    /// </summary>
    public double Lmst { get; set; }

    /// <summary>
    /// T: S_GEOM Mean anomaly, degrees
    /// </summary>
    public double MeanAnom { get; set; }

    /// <summary>
    /// T: S_GEOM Mean longitude, degrees
    /// </summary>
    public double MeanLon { get; set; }

    /// <summary>
    /// T: S_GEOM Right ascension, degrees
    /// </summary>
    public double RightAscension { get; set; }

    /// <summary>
    /// I: Surface pressure, millibars (hPa), used for refraction correction and ampress
    /// </summary>
    public double Press { get; set; }

    /// <summary>
    /// O: S_PRIME Factor that normalizes Kt, Kn, etc.
    /// </summary>
    public double Prime { get; set; }

    /// <summary>
    /// O: S_SBCF Shadow-band correction factor
    /// </summary>
    public double ShadowBandCorrFact { get; set; }

    /// <summary>
    /// I: Shadow-band width (cm)
    /// </summary>
    public double ShadowBandWidth { get; set; }

    /// <summary>
    /// I: Shadow-band radius (cm)
    /// </summary>
    public double ShadowBandRadius { get; set; }

    /// <summary>
    /// I: Shadow-band sky factor
    /// </summary>
    public double ShadowBandSkyFact { get; set; }

    /// <summary>
    /// I: Solar constant (NREL uses 1367 W/sq m)
    /// </summary>
    public double SolarConst { get; set; }

    /// <summary>
    /// T: S_SRHA Sunset(/rise) hour angle, degrees
    /// </summary>
    public double SunsetHourAng { get; set; }

    /// <summary>
    /// O: S_SRSS Sunrise time, minutes from midnight, local, WITHOUT refraction
    /// </summary>
    public double SunriseTime { get; set; }

    /// <summary>
    /// O: S_SRSS Sunset time, minutes from midnight, local, WITHOUT refraction
    /// </summary>
    public double SunsetTime { get; set; }

    /// <summary>
    /// I: Ambient dry-bulb temperature, degrees C, used for refraction correction
    /// </summary>
    public double Temp { get; set; }

    /// <summary>
    /// I: Degrees tilt from horizontal of panel
    /// </summary>
    public double PanelTilt { get; set; }

    /// <summary>
    /// I: Time zone, east (west negative). USA:  Mountain = -7, Central = -6, etc.
    /// </summary>
    public double Timezone { get; set; }

    /// <summary>
    /// T: S_TST True solar time, minutes from midnight
    /// </summary>
    public double TrueSolarTime { get; set; }

    /// <summary>
    /// T: S_TST True solar time - local standard time
    /// </summary>
    public double TrueSolarTimeLoc { get; set; }

    /// <summary>
    /// O: S_PRIME Factor that denormalizes Kt', Kn', etc.
    /// </summary>
    public double Unprime { get; set; }

    /// <summary>
    /// T: S_GEOM Universal (Greenwich) standard time
    /// </summary>
    public double TimeUtc { get; set; }

    /// <summary>
    /// T: S_ZENETR Solar zenith angle, no atmospheric correction (= ETR)
    /// </summary>
    public double SolarZenith { get; set; }

    /// <summary>
    /// O: S_REFRAC Solar zenith angle, deg. from zenith, refracted
    /// </summary>
    public double SolarZenithCorr { get; set; }
}

public readonly record struct SolarPosition(double Azimuth_deg, double Elevation_deg);
