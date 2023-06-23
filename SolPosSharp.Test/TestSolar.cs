using SolPosSharp;

namespace SolPosSharp.Test
{
    public class TestSolar
    {
        [Fact]
        public void TestSolPos()
        {
            IoData pdat = new()
            {
                Longitude = -84.43,
                Latitude = 33.65,
                Timezone = 0,
                Year = 1999,
                DayOfYear = 203, // July 22nd
                                 // UTC time
                Hour = 14,
                Minute = 45,
                Second = 37,
                Temp = 27.0,
                Press = 1006.0,
                PanelTilt = 33.65, // Tilted at latitude
                PanelAzimuth = 135.0, // 135 deg. = SE
            };
            ErrorCode retval = Solar.CalcSolarPosition(pdat); // ExSolpos function call
            TestAssertions();
            // Switch input from day of year to month/day
            pdat.Function &= ~Function.S_DOY;
            pdat.Month = 7;
            pdat.DayOfMonth = 22;
            pdat.DayOfYear = -999;
            retval = Solar.CalcSolarPosition(pdat); // ExSolpos function call
            TestAssertions();

            void TestAssertions()
            {
                // NREL    -> 1999.07.22, daynum 203, retval 0, amass 1.335752, ampress 1.326522
                Assert.Equal(1999, pdat.Year);
                Assert.Equal(7, pdat.Month);
                Assert.Equal(22, pdat.DayOfMonth);
                Assert.Equal(203, pdat.DayOfYear);
                Assert.Equal((ErrorCode)0, retval);
                Assert.Equal(1.335752, pdat.AirMass, 4);
                Assert.Equal(1.326522, pdat.AirMassPress, 4);

                // NREL    -> etr 989.668518,       etrn 1323.239868,      etrtilt 1207.547363
                Assert.Equal(989.668518, pdat.ExtTerGHI, 2);
                Assert.Equal(1323.239868, pdat.ExtTerDNI, 2);
                Assert.Equal(1207.547363, pdat.ExtTerTilt, 2);

                // NREL    -> prime 1.037040,         sbcf 1.201910,         sunrise 647.173431 (UTC)
                Assert.Equal(1.037040, pdat.Prime, 4);
                Assert.Equal(1.201910, pdat.ShadowBandCorrFact, 4);
                Assert.Equal(647.173431, pdat.SunriseTime, 2);

                // NREL    -> sunset 1481.111206 (UTC),      unprime 0.964283,          zenref 41.590069
                Assert.Equal(1481.111206, pdat.SunsetTime, 2);
                Assert.Equal(0.964283, pdat.Unprime, 5);
                Assert.Equal(41.590069, pdat.SolarZenithCorr, 3);

                Assert.Equal(97.033314, pdat.SolarAzimuth, 3);
                Assert.Equal(48.409750, pdat.SolarElevationCorr, 3);
                Assert.Equal(48.396337, pdat.SolarElevation, 3);
            }
        }
    }
}