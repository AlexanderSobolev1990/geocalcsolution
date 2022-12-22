//----------------------------------------------------------------------------------------------------------------------
///
/// \file       test_spml_geodesy.cpp
/// \brief      Тесты библиотеки spml
/// \date       01.07.21 - создан
/// \author     Соболев А.А.
///

//#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_spml_geodesy
// Boost includes:
#include <boost/test/unit_test.hpp>

// System includes:
#include <fstream>
#include <iostream>
#include <chrono>
#include <thread>

// SPML includes:
#include <geodesy.h>
//----------------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_suite_GEOtoRAD )

const SPML::Units::TAngleUnit unitAngle = SPML::Units::TAngleUnit::AU_Degree;
const SPML::Units::TRangeUnit unitRange = SPML::Units::TRangeUnit::RU_Kilometer;

const SPML::Geodesy::CEllipsoid sphere_6371 = SPML::Geodesy::Ellipsoids::Sphere6371();
const SPML::Geodesy::CEllipsoid el_WGS84 = SPML::Geodesy::Ellipsoids::WGS84();

// начальная точка (lat-lon)
double startPoint[25][2] =
{
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 }
};

// конечная точка (lat-lon)
double endPoint[25][2] =
{
    { 30.0, 0.0 },
    { 60.0, 0.0 },
    { 30.0, 30.0 },
    { 60.0, 60.0 },
    { 0.0, 30.0 },
    { 0.0, 60.0 },
    { -30.0, 30.0 },
    { -60.0, 60.0 },
    { -30.0, 0.0 },
    { -60.0, 0.0 },
    { -30.0, -30.0 },
    { -60.0, -60.0 },
    { 0.0, -30.0 },
    { 0.0, -60.0 },
    { 30.0, -30.0 },
    { 60.0, -60.0 },
    { 59.9386300, 30.3141300 },     // Saint-Petersburg, Russia
    { 52.5243700, 13.4105300 },     // Berlin, Germany
    { 51.5085300, -0.1257400 },     // London, Great Britain
    { 48.8534100, 2.3488000 },      // Paris, France
    { 41.8919300, 12.5113300 },     // Roma, Italy
    { 64.1354800, -21.8954100 },    // Reykjavik, Iceland
    { 40.7142700, -74.0059700 },    // New-York, USA
    { -34.6131500, -58.3772300 },   // Buenos Aires, Argentina
    { -33.8678500, 151.2073200 }    // Sydney, Australia
};

double rightAnswerInverseProblem_Sphere6371000[25][3] =
{
    { 3335.847799336762, 0.0, 180.0 },
    { 6671.695598673524, 0.0, 180.0 },
    { 4604.53989281927, 40.893394649131, 229.106605350869 },
    { 8397.717492500104, 26.565051177078, 243.434948822922 },
    { 3335.847799336762, 90.0, 270.0 },
    { 6671.695598673524, 90.0, 270.0 },
    { 4604.53989281927, 139.106605350869, 310.893394649131 },
    { 8397.717492500104, 153.434948822922, 296.565051177078 },
    { 3335.847799336762, 180.0, 0.0 },
    { 6671.695598673524, 180.0, 0.0 },
    { 4604.53989281927, 220.893394649131, 49.106605350869 },
    { 8397.717492500104, 206.565051177078, 63.434948822922 },
    { 3335.847799336762, 270.0, 90.0 },
    { 6671.695598673524, 270.0, 90.0 },
    { 4604.53989281927, 319.106605350869, 130.893394649131 },
    { 8397.717492500104, 333.434948822922, 116.565051177078 },

    { 634.430920264603, 320.181220542521, 133.993227500364 },
    { 1608.171539213493, 267.225012356327, 67.500629141408 },
    { 2500.279541412505, 275.046184204728, 64.249739712143 },
    { 2486.515274556244, 266.942048390901, 58.657801947884 },
    { 2375.139101545627, 240.124400843247, 40.960381828756 },
    { 3305.826756198104, 310.708086770129, 77.933216738675 },
    { 7510.302978483882, 310.317895818262, 34.47939646295 },
    { 13475.866954006504, 253.102223319853, 40.864960821945 },
    { 14496.045468976658, 92.928775544395, 317.398997046878 }
};

double rightAnswerDirectProblem_Sphere6371000[25][3] =
{
    { 30.0, 0.0, 0.0 },
    { 60.0, 0.0, 0.0 },
    { 30.00000, 30.00000, 49.10661 },
    { 60.00000, 60.00000, 63.43495 },
    { 0.00000, 30.00000, 90.00000 },
    { 0.00000, 60.00000, 90.00000 },
    { -30.00000, 30.00000, 130.89339 },
    { -60.00000, 60.00000, 116.56505 },
    { -30.00000, 0.00000, 180.0 },
    { -60.00000, 0.00000, 180.0 },
    { -30.00000, -30.00000, 229.10661 },
    { -60.00000, -60.00000, 243.43495 },
    { 0.00000, -30.00000, 270.00000 },
    { 0.00000, -60.00000, 270.00000 },
    { 30.00000, -30.00000, 310.89339 },
    { 60.00000, -60.00000, 296.56505 }, // 16

    { 59.93863, 30.31413, 313.99323 },
    { 52.52437, 13.41053, 247.50063 },
    { 51.50853, -0.12574, 244.24974 },
    { 48.85341, 2.34880, 238.65780 }, // 20

    { 41.89193, 12.51133, 220.96038 },
    { 64.13548, -21.89541, 257.93322 },
    { 40.71427, -74.00597, 214.47940 },
    { -34.61315, -58.37723, 220.86496 },
    { -33.86785, 151.20732, 137.39900 } // 25
};

double rightAnswerInverseProblem_WGS84[25][3] =
{
    { 3320.11338490278, 0.0, 180.0 },
    { 6654.07283255148, 0.0, 180.0 },
    { 4596.22309715941, 41.066728349773, 229.282146152291 },
    { 8389.65383524219, 26.6683364994786, 243.559103061872 },
    { 3339.5847237982, 90.0, 270.0 },
    { 6679.16944759641, 90.0, 270.0 },
    { 4596.22309715941, 138.933271650227, 310.717853847709 },
    { 8389.65383524219, 153.331663500521, 296.440896938128 },
    { 3320.11338490278, 180.0, 0.0 },
    { 6654.07283255148, 180.0, 0.0 },
    { 4596.22309715941, 221.066728349773, 49.2821461522913 },
    { 8389.65383524219, 206.668336499479, 63.5591030618716 },
    { 3339.5847237982, 270.0, 90.0 },
    { 6679.16944759641, 270.0, 90.0 },
    { 4596.22309715941, 318.933271650227, 130.717853847709 },
    { 8389.65383524219, 333.331663500521, 116.440896938128 }, // 16

    { 636.015930287104, 320.126944770617, 133.93894388736 },
    { 1613.33708756292, 267.253437457649, 67.5288169619985 },
    { 2508.31301483815, 275.070800591577, 64.2734524345967 },
    { 2493.92247471802, 266.983302187342, 58.698181706953 }, // 20

    { 2379.36533650612, 240.204880398643, 41.0402378345664 },
    { 3317.37231308908, 310.685756113135, 77.9091383324392 },
    { 7531.172180035, 310.353662800418, 34.492624273887 },
    { 13459.2791448219, 253.278786758906, 40.9713451154643 },
    { 14484.0079887903, 92.7281133803279, 317.323828722842 } //25
};

double rightAnswerDirectProblem_WGS84[25][3] =
{
    { 30.0, 0.0, 0.0 },
    { 60.0, 0.0, 0.0 },
    { 30.00000, 30.00000, 49.282146445946 },
    { 60.00000, 60.00000, 63.5590965976542 },
    { 0.00000, 30.00000, 90.00000 },
    { 0.00000, 60.00000, 90.00000 },
    { -30.00000, 30.00000, 130.717853554054 },
    { -60.00000, 60.00000, 116.440905302634 },
    { -30.00000, 0.00000, 180.0 },
    { -60.00000, 0.00000, 180.0 },
    { -30.00000, -30.00000, 229.282146445946 },
    { -60.00000, -60.00000, 243.559094697366 },
    { 0.00000, -30.00000, 270.00000 },
    { 0.00000, -60.00000, 270.00000 },
    { 30.00000, -30.00000, 310.717838319556 },
    { 60.00000, -60.00000, 296.440905302634 }, // 16

    { 59.9386300, 30.3141300, 313.938959560895 },
    { 52.5243700, 13.4105300, 247.528820831467 },
    { 51.5085300, -0.1257400, 244.273443839833 },
    { 48.8534100, 2.3488000, 238.698183934423 }, // 20

    { 41.8919300, 12.5113300, 221.040246161089 },
    { 64.1354800, -21.8954100, 257.9091202346 },
    { 40.7142700, -74.0059700, 214.492623699066 },
    { -34.6131500, -58.3772300, 220.971345375481 },
    { -33.8678500, 151.2073200, 137.323827838602 } // 25
};

BOOST_AUTO_TEST_CASE( test_GEOtoRAD_1 )
{
    const int i = 0;
    double eps = 0.00001;
    SPML::Geodesy::RAD rad = SPML::Geodesy::GEOtoRAD( el_WGS84, unitRange, unitAngle,
        SPML::Geodesy::Geographic( startPoint[i][0], startPoint[i][1] ),
        SPML::Geodesy::Geographic( endPoint[i][0], endPoint[i][1] ) );
    BOOST_CHECK_CLOSE( rad.R, rightAnswerInverseProblem_WGS84[i][0], eps );
    BOOST_CHECK_CLOSE( rad.Az, rightAnswerInverseProblem_WGS84[i][1], eps );
    BOOST_CHECK_CLOSE( rad.AzEnd, SPML::Convert::AngleTo360( rightAnswerInverseProblem_WGS84[i][2] + 180,
        SPML::Units::AU_Degree ), eps );
}

BOOST_AUTO_TEST_CASE( test_GEOtoRAD_2 )
{
    const int i = 0;
    double eps = 0.00001;
    double R, Az, AzEnd;
    SPML::Geodesy::GEOtoRAD( el_WGS84, unitRange, unitAngle,
        startPoint[i][0], startPoint[i][1], endPoint[i][0], endPoint[i][1], R, Az, AzEnd );
    BOOST_CHECK_CLOSE( R, rightAnswerInverseProblem_WGS84[i][0], eps );
    BOOST_CHECK_CLOSE( Az, rightAnswerInverseProblem_WGS84[i][1], eps );
    BOOST_CHECK_CLOSE( AzEnd, SPML::Convert::AngleTo360( rightAnswerInverseProblem_WGS84[i][2] + 180,
        SPML::Units::AU_Degree ), eps );
}

BOOST_AUTO_TEST_CASE( test_GEOtoRAD_RADtoGEO_Sphere )
{
    std::ofstream out("test_GEOtoRAD_RADtoGEO.txt");

    double latEnd;
    double lonEnd;
    double range;
    double az;
    double az2;
    double eps = 0.0001;
    double small = 1.0e-1;

    for( int i = 0; i < 25; i++ ) {
        SPML::Geodesy::GEOtoRAD( sphere_6371, unitRange, unitAngle,
            startPoint[i][0], startPoint[i][1],
            endPoint[i][0], endPoint[i][1],
            range, az, az2 );
        out << "----------------------------------------------------------------------------------------------" << std::endl;
        out << "Test #" << ( i + 1 ) << std::endl;
        out << "StartPoint:" << std::endl;
        out << startPoint[i][0] << ' ' << startPoint[i][1] << std::endl;
        out << "EndPoint:" << std::endl;
        out << endPoint[i][0] << ' ' << endPoint[i][1] << std::endl;

        out << "Result Inverse Problem:" << std::endl;
        out << range << ' ' << az << ' ' << az2 << std::endl;
        BOOST_CHECK_CLOSE_FRACTION( range, rightAnswerInverseProblem_Sphere6371000[i][0], eps );
        BOOST_CHECK_CLOSE_FRACTION( az, rightAnswerInverseProblem_Sphere6371000[i][1], eps );
        BOOST_CHECK_CLOSE_FRACTION( az2, SPML::Convert::AngleTo360( rightAnswerInverseProblem_Sphere6371000[i][2] + 180.0,
            SPML::Units::TAngleUnit::AU_Degree ), eps );

        SPML::Geodesy::RADtoGEO( sphere_6371, unitRange, unitAngle,
            startPoint[i][0], startPoint[i][1],
            range, az,
            latEnd, lonEnd, az2 );

        out << "Result Direct Problem:" << std::endl;
        out << latEnd << ' ' << lonEnd << ' ' << az2 << std::endl;
//        out << "RightAnswer Direct Problem:" << std::endl;
//        out << rightAnswerDirectProblem_Sphere6371000[i][0] << ' ' << rightAnswerDirectProblem_Sphere6371000[i][1] << ' ' << rightAnswerDirectProblem_Sphere6371000[i][2] << std::endl;
        BOOST_CHECK( std::abs( latEnd - rightAnswerDirectProblem_Sphere6371000[i][0] ) < eps );
        BOOST_CHECK( std::abs( lonEnd - rightAnswerDirectProblem_Sphere6371000[i][1] ) < eps );
        BOOST_CHECK( std::abs( az2 - rightAnswerDirectProblem_Sphere6371000[i][2] ) < eps );
    }
    out.close();
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
const SPML::Units::TAngleUnit unitAngle = SPML::Units::TAngleUnit::AU_Degree;
const SPML::Units::TRangeUnit unitRange = SPML::Units::TRangeUnit::RU_Meter;
const SPML::Geodesy::CEllipsoid el = SPML::Geodesy::Ellipsoids::WGS84();
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_GEOtoECEF )

double lat = 57.02929569; // град
double lon = 9.950248114; // град
double h = 56.95; // м
double x = 3426949.397; // м
double y = 601195.852; // м
double z = 5327723.994; // м

const double eps = 1.0e-3;

BOOST_AUTO_TEST_CASE( test_GEOtoECEF )
{
    double x_, y_, z_;
    SPML::Geodesy::GEOtoECEF( el, unitRange, unitAngle, lat, lon, h, x_, y_, z_ );
    BOOST_CHECK_CLOSE_FRACTION( x, x_, eps );
    BOOST_CHECK_CLOSE_FRACTION( y, y_, eps );
    BOOST_CHECK_CLOSE_FRACTION( z, z_, eps );
}

BOOST_AUTO_TEST_CASE( test_ECEFtoGEO )
{
    double lat_, lon_, h_;
//    SPML::Geodesy::ECEFtoGEO( el, unitRange, unitAngle, x, y, z, lat_, lon_, h_ );

    SPML::Geodesy::latlon( x, y, z, &lat_, &lon_, &h_ );

    BOOST_CHECK_CLOSE_FRACTION( lat, lat_, eps );
    BOOST_CHECK_CLOSE_FRACTION( lon, lon_, eps );
    BOOST_CHECK_CLOSE_FRACTION( h, h_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_GEOtoENU )

double lat = 42.002582; // град
double lon = -81.997752; // град
double h = 1139.7; // м
double lat0 = 42.0; // град
double lon0 = -82.0; // град
double h0 = 200.0; // м
double e = 186.28; // м
double n = 286.84; // м
double u = 939.69; // м

const double eps = 1.0e-3;

BOOST_AUTO_TEST_CASE( test_GEOtoENU )
{
    double e_, n_, u_;
    SPML::Geodesy::GEOtoENU( el, unitRange, unitAngle, lat, lon, h, lat0, lon0, h0, e_, n_, u_ );
    BOOST_CHECK_CLOSE_FRACTION( e, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( n, n_, eps );
    BOOST_CHECK_CLOSE_FRACTION( u, u_, eps );
}

BOOST_AUTO_TEST_CASE( test_ENUtoGEO )
{
    double lat_, lon_, h_;
    SPML::Geodesy::ENUtoGEO( el, unitRange, unitAngle, e, n, u, lat0, lon0, h0, lat_, lon_, h_ );
    BOOST_CHECK_CLOSE_FRACTION( lat, lat_, eps );
    BOOST_CHECK_CLOSE_FRACTION( lon, lon_, eps );
    BOOST_CHECK_CLOSE_FRACTION( h, h_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_GEOtoAER )

double lat = 42.002582; // град
double lon = -81.997752; // град
double h = 1139.7; // м
double lat0 = 42.0; // град
double lon0 = -82.0; // град
double h0 = 200.0; // м
double a = 33.0; // м
double e = 70.0; // м
double r = 1000.0; // м

const double eps = 1.0e-3;

BOOST_AUTO_TEST_CASE( test_GEOtoAER )
{
    double a_, e_, r_;
    SPML::Geodesy::GEOtoAER( el, unitRange, unitAngle, lat, lon, h, lat0, lon0, h0, a_, e_, r_ );
    BOOST_CHECK_CLOSE_FRACTION( a, a_, eps );
    BOOST_CHECK_CLOSE_FRACTION( e, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( r, r_, eps );
}

BOOST_AUTO_TEST_CASE( test_AERtoGEO )
{
    double lat_, lon_, h_;
    SPML::Geodesy::AERtoGEO( el, unitRange, unitAngle, a, e, r, lat0, lon0, h0, lat_, lon_, h_ );
    BOOST_CHECK_CLOSE_FRACTION( lat, lat_, eps );
    BOOST_CHECK_CLOSE_FRACTION( lon, lon_, eps );
    BOOST_CHECK_CLOSE_FRACTION( h, h_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_ECEFtoENU )

double x = 660930.192761082;
double y = -4701424.222957011;
double z = 4246579.604632881;
double lat0 = 42.0; // град
double lon0 = -82.0; // град
double h0 = 200.0; // м
double e = 186.27752; // м
double n = 286.84222; // м
double u = 939.69262; // м

double eps = 1.0e-5;

BOOST_AUTO_TEST_CASE( test_ECEFtoENU )
{
    double e_, n_, u_;
    SPML::Geodesy::ECEFtoENU( el, unitRange, unitAngle, x, y, z, lat0, lon0, h0, e_, n_, u_ );
    BOOST_CHECK_CLOSE_FRACTION( e, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( n, n_, eps );
    BOOST_CHECK_CLOSE_FRACTION( u, u_, eps );
}

BOOST_AUTO_TEST_CASE( test_ENUtoECEF )
{
    double x_, y_, z_;
    SPML::Geodesy::ENUtoECEF( el, unitRange, unitAngle, e, n, u, lat0, lon0, h0, x_, y_, z_ );
    BOOST_CHECK_CLOSE_FRACTION( x, x_, eps );
    BOOST_CHECK_CLOSE_FRACTION( y, y_, eps );
    BOOST_CHECK_CLOSE_FRACTION( z, z_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_ECEFtoAER )

double x = 660930.192761082;
double y = -4701424.222957011;
double z = 4246579.604632881;
double lat0 = 42.0; // град
double lon0 = -82.0; // град
double h0 = 200.0; // м
double a = 33.0; // м
double e = 70.0; // м
double r = 1000.0; // м

double eps = 1.0e-5;

BOOST_AUTO_TEST_CASE( test_ECEFtoAER )
{
    double a_, e_, r_;
    SPML::Geodesy::ECEFtoAER( el, unitRange, unitAngle, x, y, z, lat0, lon0, h0, a_, e_, r_ );
    BOOST_CHECK_CLOSE_FRACTION( a, a_, eps );
    BOOST_CHECK_CLOSE_FRACTION( e, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( r, r_, eps );
}

BOOST_AUTO_TEST_CASE( test_AERtoECEF )
{
    double x_, y_, z_;
    SPML::Geodesy::AERtoECEF( el, unitRange, unitAngle, a, e, r, lat0, lon0, h0, x_, y_, z_ );
    BOOST_CHECK_CLOSE_FRACTION( x, x_, eps );
    BOOST_CHECK_CLOSE_FRACTION( y, y_, eps );
    BOOST_CHECK_CLOSE_FRACTION( z, z_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE( test_suite_ENUtoAER )

double e = 186.27752; // м
double n = 286.84222; // м
double u = 939.69262; // м
double a1 = 33.0; // м
double e1 = 70.0; // м
double r1 = 1000.0; // м

double eps = 1.0e-5;

BOOST_AUTO_TEST_CASE( test_ENUtoAER )
{
    double a_, e_, r_;
    SPML::Geodesy::ENUtoAER( unitRange, unitAngle, e, n, u, a_, e_, r_ );
    BOOST_CHECK_CLOSE_FRACTION( a1, a_, eps );
    BOOST_CHECK_CLOSE_FRACTION( e1, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( r1, r_, eps );
}

BOOST_AUTO_TEST_CASE( test_AERtoENU )
{
    double e_, n_, u_;
    SPML::Geodesy::AERtoENU( unitRange, unitAngle, a1, e1, r1, e_, n_, u_ );
    BOOST_CHECK_CLOSE_FRACTION( e, e_, eps );
    BOOST_CHECK_CLOSE_FRACTION( n, n_, eps );
    BOOST_CHECK_CLOSE_FRACTION( u, u_, eps );
}

BOOST_AUTO_TEST_SUITE_END()
//----------------------------------------------------------------------------------------------------------------------


