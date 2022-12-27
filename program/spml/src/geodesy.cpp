//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.cpp
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <geodesy.h>

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{
namespace Geodesy /// Геодезические функции и функции перевода координат
{

//----------------------------------------------------------------------------------------------------------------------
CEllipsoid::CEllipsoid()
{
    a = 0.0;
    b = 0.0;
    f = 0.0;
    invf = 0.0;
}

CEllipsoid::CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef )
{
    name = ellipsoidName;
    a = semiMajorAxis;
    invf = inverseFlattening;
    if( isInvfDef && ( Compare::IsZeroAbs( inverseFlattening ) || std::isinf( inverseFlattening ) ) ) {
        b = semiMajorAxis;
        f = 0.0;
    } else if ( isInvfDef ) {
        b = ( 1.0 - ( 1.0 / inverseFlattening ) ) * semiMajorAxis;
        f = 1.0 / inverseFlattening;
    } else {
        b = semiMinorAxis;
        f = 1.0 / inverseFlattening;
    }
}

//CEllipsoid::CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening )
//{
//    name = ellipsoidName;
//    a = semiMajorAxis;
//    invf = inverseFlattening;
//    if( Compare::IsZero( semiMinorAxis ) && !Compare::IsZero( inverseFlattening ) ) { // b = 0, invf != 0
//        b = ( 1.0 - ( 1.0 / inverseFlattening ) ) * semiMajorAxis;
//        f = 1.0 / inverseFlattening;
//    } else if( !Compare::IsZero( semiMinorAxis ) && Compare::IsZero( inverseFlattening ) ) {
//        b = semiMinorAxis;
//        f = 0.0;
//    } else {
//        b = semiMinorAxis;
//        f = 1.0 / inverseFlattening;
//    }
//}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // По умолчанию Радианы:
    double _latStart = latStart;
    double _lonStart = lonStart;
    double _latEnd = latEnd;
    double _lonEnd = lonEnd;

    // При необходимости переведем входные данные в Радианы:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _latEnd *= Convert::DgToRdD;
            _lonEnd *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqualAbs( a, b ) ) { // При расчете на сфере используем упрощенные формулы
        // Azimuth
        double fact1, fact2, fact3;
        fact1 = std::cos( _latEnd ) * std::sin( _lonEnd - _lonStart );
        fact2 = std::cos( _latStart ) * std::sin( _latEnd );
        fact3 = std::sin( _latStart ) * std::cos( _latEnd ) * std::cos( _lonEnd - _lonStart );
        az = Convert::AngleTo360( std::atan2( fact1, fact2 - fact3 ), Units::AU_Radian ); // [рад] - Прямой азимут в начальной точке

        // ReverseAzimuth
        fact1 = std::cos( _latStart ) * std::sin( _lonEnd - _lonStart );
        fact2 = std::cos( _latStart ) * std::sin( _latEnd ) * std::cos( _lonEnd - _lonStart );
        fact3 = std::sin( _latStart ) * std::cos( _latEnd );
        azEnd = Convert::AngleTo360( ( std::atan2( fact1, fact2 - fact3 ) ), Units::AU_Radian ); // [рад] - Прямой азимут в конечной точке

        // Distance
        double temp1, temp2, temp3;
        temp1 = std::sin( _latStart ) * std::sin( _latEnd );
        temp2 = std::cos( _latStart ) * std::cos( _latEnd ) * std::cos( _lonEnd - _lonStart );
        temp3 = temp1 + temp2;        
        d = std::acos( temp3 ) * a ; // [м]
    } else { // Для эллипсоида используем формулы Винсента
        double L = _lonEnd - _lonStart;

        double U1 = std::atan( ( 1.0 - f ) * std::tan( _latStart ) );
        double U2 = std::atan( ( 1.0 - f ) * std::tan( _latEnd ) );

        double sinU1 = std::sin( U1 );
        double cosU1 = std::cos( U1 );
        double sinU2 = std::sin( U2 );
        double cosU2 = std::cos( U2 );

        // eq. 13
        double lambda = L;
        double lambda_new = 0.0;
        int iterLimit = 100;

        double sinSigma = 0.0;
        double cosSigma = 0.0;
        double sigma = 0.0;
        double sinAlpha = 0.0;
        double cosSqAlpha = 0.0;
        double cos2SigmaM = 0.0;
        double c = 0.0;
        double sinLambda = 0.0;
        double cosLambda = 0.0;

        do {
            sinLambda = std::sin( lambda );
            cosLambda = std::cos( lambda );

            // eq. 14
            sinSigma = std::sqrt( ( ( cosU2 * sinLambda ) * ( cosU2 * sinLambda ) +
                ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) * ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ) );
            if( Compare::IsZeroAbs( sinSigma ) ) { // co-incident points
                d = 0.0;
                az = 0.0;
                azEnd = 0.0;
                return;
            }

            // eq. 15
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

            // eq. 16
            sigma = std::atan2( sinSigma, cosSigma );

            // eq. 17    Careful!  sin2sigma might be almost 0!
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            cosSqAlpha = 1 - sinAlpha * sinAlpha;

            // eq. 18    Careful!  cos2alpha might be almost 0!
            cos2SigmaM = cosSigma - 2.0 * sinU1 * sinU2 / cosSqAlpha;

            if( std::isnan( cos2SigmaM ) ) {
                cos2SigmaM = 0; // equatorial line: cosSqAlpha = 0
            }

            // eq. 10
            c = ( f / 16.0 ) * cosSqAlpha * ( 4.0 + f * ( 4.0 - 3.0 * cosSqAlpha ) );

            lambda_new = lambda;

            // eq. 11 (modified)
            lambda = L + ( 1.0 - c ) * f * sinAlpha *
                ( sigma + c * sinSigma * ( cos2SigmaM + c * cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) ) );

        } while( std::abs( ( lambda - lambda_new ) / lambda ) > 1.0e-15 && --iterLimit > 0 ); // see how much improvement we got

        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + uSq / 16384.0 * ( 4096.0 + uSq * ( -768.0 + uSq * ( 320.0 - 175.0 * uSq ) ) );

        // eq. 4
        double B = uSq / 1024.0 * ( 256.0 + uSq * ( -128.0 + uSq * ( 74.0 - 47.0 * uSq ) ) );

        // eq. 6
        double deltaSigma = B * sinSigma *
            ( cos2SigmaM + ( B / 4.0 ) * ( cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) -
            ( B / 6.0 ) * cos2SigmaM * ( -3.0 + 4.0 * sinSigma * sinSigma ) * ( -3.0 + 4.0 * cos2SigmaM * cos2SigmaM ) ) );

        // eq. 19        
        d = b * A * ( sigma - deltaSigma ); // [m]

        // eq. 20
        az = Convert::AngleTo360( std::atan2( ( cosU2 * sinLambda ),
            ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в начальной точке, [рад]

        // eq. 21
        azEnd = Convert::AngleTo360( std::atan2( ( cosU1 * sinLambda ), ( -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // az, azEnd, d сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgD;
            azEnd *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            d *= 0.001;
            break;
        }
        default:
            assert( false );
    }
    return;
}

RAD GEOtoRAD(const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const Geographic &end )
{
    double d, az, azEnd;
    GEOtoRAD( ellipsoid, rangeUnit, angleUnit, start.Lat, start.Lon, end.Lat, end.Lon, d, az, azEnd );
    return RAD( d, az, azEnd );
}

//----------------------------------------------------------------------------------------------------------------------
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // по умолчанию Метры-Радианы:
    double _latStart = latStart;    // [рад]
    double _lonStart = lonStart;    // [рад]
    double _d = d;                  // [м]
    double _az = az;                // [рад]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _az *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _d *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqualAbs(a, b) ) { // При расчете на сфере используем упрощенные формулы
        _d = _d / a; // Нормирование
        // latitude
        double temp1, temp2, temp3;
        temp1 = std::sin( _latStart ) * std::cos( _d );
        temp2 = std::cos( _latStart ) * std::sin( _d ) * std::cos( _az );
        latEnd = std::asin( temp1 + temp2 ); // [рад]

        // longitude
        temp1 = std::sin( _d ) * std::sin( _az );
        temp2 = std::cos( _latStart ) * std::cos( _d );
        temp3 = std::sin( _latStart ) * std::sin( _d ) * std::cos( _az );
        lonEnd = _lonStart + std::atan2( temp1, temp2 - temp3 ); // [рад]

        // final bearing
        temp1 = std::cos( _latStart ) * std::sin( _az );
        temp2 = std::cos( _latStart ) * std::cos( _d ) * std::cos( _az );
        temp3 = std::sin( _latStart ) * std::sin( _d );
        azEnd = Convert::AngleTo360( std::atan2( temp1, temp2 - temp3 ), Units::TAngleUnit::AU_Radian ); // [рад] - Прямой азимут в конечной точке
    } else { // Для эллипсоида используем формулы Винсента
        double cosAlpha1 = std::cos( _az );
        double sinAlpha1 = std::sin( _az );
        double s = _d; // distance [m]
        double tanU1 = ( 1.0 - f ) * std::tan( _latStart );
        double cosU1 = 1.0 / std::sqrt( ( 1.0 + tanU1 * tanU1 ) );
        double sinU1 = tanU1 * cosU1;

        // eq. 1
        double sigma1 = std::atan2( tanU1, cosAlpha1 );

        // eq. 2
        double sinAlpha = cosU1 * sinAlpha1;
        double cosSqAlpha = 1 - sinAlpha * sinAlpha;
        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1.0 + ( uSq / 16384.0 ) * ( 4096.0 + uSq * ( -768.0 + uSq * ( 320.0 - 175.0 * uSq ) ) );

        // eq. 4
        double B = ( uSq / 1024.0 ) * ( 256.0 + uSq * ( -128.0 + uSq * ( 74.0 - 47.0 * uSq ) ) );

        // iterate until there is a negligible change in sigma
        double sOverbA = s / ( b * A );
        double sigma = sOverbA;
        double prevSigma = sOverbA;
        double cos2SigmaM = 0.0;
        double sinSigma = 0.0;
        double cosSigma = 0.0;
        double deltaSigma = 0.0;

        int iterations = 0;

        while( true ) {
            // eq. 5
            cos2SigmaM = std::cos( 2.0 * sigma1 + sigma );
            sinSigma = std::sin( sigma );
            cosSigma = std::cos( sigma );

            // eq. 6
            deltaSigma = B * sinSigma * ( cos2SigmaM +
                ( B / 4.0 ) * ( cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) -
                ( B / 6.0 ) * cos2SigmaM * ( -3.0 + 4.0 * sinSigma * sinSigma ) * ( -3.0 + 4.0 * cos2SigmaM * cos2SigmaM ) ) );

            // eq. 7
            sigma = sOverbA + deltaSigma;

            // break after converging to tolerance
            if( std::abs( sigma - prevSigma ) < 1.0e-15 || std::isnan( std::abs( sigma - prevSigma ) ) ) {
                break;
            }
            prevSigma = sigma;

            iterations++;
            if( iterations > 1000 ) {
                break;
            }
        }
        cos2SigmaM = std::cos( 2.0 * sigma1 + sigma );
        sinSigma = std::sin( sigma );
        cosSigma = std::cos( sigma );

        double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;

        // eq. 8
        latEnd = std::atan2( sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
            ( 1.0 - f ) * std::sqrt( ( sinAlpha * sinAlpha + tmp * tmp) ) ); // [рад]

        // eq. 9
        double lambda = std::atan2( ( sinSigma * sinAlpha1 ), ( cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1 ) );

        // eq. 10
        double c = ( f / 16.0 ) * cosSqAlpha * ( 4.0 + f * ( 4.0 - 3.0 * cosSqAlpha ) );

        // eq. 11
        double L = lambda - ( 1.0 - c ) * f * sinAlpha * ( sigma + c * sinSigma *
            ( cos2SigmaM + c * cosSigma * ( -1.0 + 2.0 * cos2SigmaM * cos2SigmaM ) ) );

        //double phi = ( _lonStart + L + 3 * PI ) % ( 2 * PI ) - PI;   // to -180.. 180 original!
        //lonEnd = ( _lonStart + L ) * RdToDgD;   // [град] my
        lonEnd = _lonStart + L;   // [рад] my

        // eq. 12
        double alpha2 = std::atan2( sinAlpha, -tmp ); // final bearing, if required
        azEnd = Convert::AngleTo360( alpha2, Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // latEnd, lonEnd, azEnd сейчас в радианах

    // Проверим, нужен ли перевод:    
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            latEnd *= Convert::RdToDgD;
            lonEnd *= Convert::RdToDgD;
            azEnd *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    return;
}

Geographic RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const RAD &rad, double &azEnd )
{
    double latEnd, lonEnd;
    RADtoGEO( ellipsoid, rangeUnit, angleUnit,
        start.Lat, start.Lon, rad.R, rad.Az, latEnd, lonEnd, azEnd );
    return Geographic( latEnd, lonEnd );
}
//----------------------------------------------------------------------------------------------------------------------
void GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double &x, double &y, double &z )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
//    double b = ellipsoid.B();

    // по умолчанию Метры-Радианы:
    double _lat = lat; // [рад]
    double _lon = lon; // [рад]
    double _h = h;     // [м]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRdD;
            _lon *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце
    assert( !Compare::IsZeroAbs( a * a ) );
//    double es = 1.0 - ( ( b * b ) / ( a * a ) ); // e^2
    double es = ellipsoid.EccentricityFirstSquared();

    double sinLat = std::sin( _lat );
    double cosLat = std::cos( _lat );
    double sinLon = std::sin( _lon );
    double cosLon = std::cos( _lon );

    double arg = 1.0 - ( es * ( sinLat * sinLat ) );
    assert( arg > 0 );
    double v = a / std::sqrt( arg );

    x = ( v + _h ) * cosLat * cosLon; // [м]
    y = ( v + _h ) * cosLat * sinLon; // [м]
    z = ( v * ( 1.0 - es ) + _h ) * std::sin( _lat );    // [м]

    // x, y, x сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            x /= 1000.0;
            y /= 1000.0;
            z /= 1000.0;
            break;
        }
        default:
            assert( false );
    }
}

XYZ GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic point )
{
    double x, y, z;
    GEOtoECEF( ellipsoid, rangeUnit, angleUnit, point.Lat, point.Lon, point.Height, x, y, z );
    return XYZ( x, y, z );
}
//----------------------------------------------------------------------------------------------------------------------
/*
// OLD METHOD
void ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double &lat, double &lon, double &h )
{
    // An Improved Algorithm for Geocentric to Geodetic Coordinate Conversion, Ralph Toms, Feb 1996. UCRL-JC-123138

    static const double AD_C = 1.0026000;   // Toms region 1 constant. ( h_min = -1e5 [m], h_max = 2e6 [m] ) - MOST CASES!
//    static const double AD_C = 1.00092592;  // Toms region 2 constant. ( h_min = 2e6 [m], h_max = 6e6 [m] )
//    static const double AD_C = 0.999250297; // Toms region 3 constant. ( h_min = 6e6 [m], h_max = 18e6 [m] )
//    static const double AD_C = 0.997523508; // Toms region 4 constant. ( h_min = 18e5 [m], h_max = 1e9 [m] )
    static const double COS_67P5 = 0.38268343236508977; // Cosine of 67.5 degrees

    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();

//    double es = 1.0 - ( b * b ) / ( a * a );    // Eccentricity squared : (a^2 - b^2)/a^2
//    double ses = ( a * a ) / ( b * b ) - 1.0;   // Second eccentricity squared : (a^2 - b^2)/b^2
    double es = ellipsoid.EccentricityFirstSquared();   // Eccentricity squared : (a^2 - b^2)/a^2
    double ses = ellipsoid.EccentricitySecondSquared(); // Second eccentricity squared : (a^2 - b^2)/b^2

    bool At_Pole = false; // indicates whether location is in polar region

    double _x = x; // [м]
    double _y = y; // [м]
    double _z = z; // [м]

    // При необходимости переведем в Метры:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _x *= 1000.0;
            _y *= 1000.0;
            _z *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце
    double _lon = 0;
    double _lat = 0;
    double _h = 0;

    if( !Compare::IsZeroAbs( x ) ) { //if (x != 0.0)
        _lon = std::atan2( _y, _x );
    } else {
        if( _y > 0 ) {
            _lon = Consts::PI_D / 2;
        } else if ( _y < 0 ) {
            _lon = -Consts::PI_D * 0.5;
        } else {
            At_Pole = true;
            _lon = 0.0;
            if( _z > 0.0 ) {            // north pole
                _lat = Consts::PI_D * 0.5;
            } else if ( _z < 0.0 ) {    // south pole
                _lat = -Consts::PI_D * 0.5;
            } else {                    // center of earth
                //lat = PI_D * 0.5 * RdToDgD; // [град]
                lat = Consts::PI_D * 0.5; // [рад] TODO: Как тут улучшить?
                lon = 0;
                h = -b;
                return;
            }
        }
    }
    double W2 = _x * _x + _y * _y; // Square of distance from Z axis
    assert( W2 > 0 );
    double W = std::sqrt( W2 ); // distance from Z axis
    double T0 = _z * AD_C; // initial estimate of vertical component
    assert( ( T0 * T0 + W2 ) > 0 );
    double S0 = std::sqrt( T0 * T0 + W2 ); //initial estimate of horizontal component
    double Sin_B0 = T0 / S0; //std::sin(B0), B0 is estimate of Bowring aux variable
    double Cos_B0 = W / S0; //std::cos(B0)
    double Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;//Math.Pow(Sin_B0, 3);
    double T1 = _z + b * ses * Sin3_B0; //corrected estimate of vertical component
    double Sum = W - a * es * Cos_B0 * Cos_B0 * Cos_B0; //numerator of std::cos(phi1)
    assert( ( T1 * T1 + Sum * Sum ) > 0 );
    double S1 = std::sqrt( T1 * T1 + Sum * Sum ); //corrected estimate of horizontal component
    double Sin_p1 = T1 / S1; //std::sin(phi1), phi1 is estimated latitude
    double Cos_p1 = Sum / S1; //std::cos(phi1)
    assert( ( 1.0 - es * Sin_p1 * Sin_p1 ) > 0 );
    double Rn = a / std::sqrt( 1.0 - es * Sin_p1 * Sin_p1 ); //Earth radius at location
    if( Cos_p1 >= COS_67P5 ) {
        _h = W / Cos_p1 - Rn;
    } else if ( Cos_p1 <= -COS_67P5 ) {
        assert( !Compare::IsZeroAbs( Cos_p1 ) );
        _h = W / ( -Cos_p1 ) - Rn;
    } else {
        assert( !Compare::IsZeroAbs( Sin_p1 ) );
        _h = ( _z / Sin_p1 ) + Rn * ( es - 1.0 );
    }
    if( !At_Pole ) {
        assert( !Compare::IsZeroAbs( Cos_p1 ) );
        _lat = std::atan2( Sin_p1, Cos_p1 );
    }
    lat = _lat; // [рад]
    lon = _lon; // [рад]
    h = _h;     // [м]

    // LLH сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat *= Convert::RdToDgD;
            lon *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            h /= 1000.0;
            break;
        }
        default:
            assert( false );
    }
}
*/

void ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double &lat, double &lon, double &h )
{
    // Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476. https://doi.org/10.1109/7.481290

    // Параметры эллипсоида:
    double _a = ellipsoid.A();
    double _es = ellipsoid.EccentricityFirstSquared();   // Eccentricity squared : (a^2 - b^2)/a^2
    const double _a1 = _a * _es;
    const double _a2 = _a1 * _a1;
    const double _a3 = _a1 * _es / 2.0;
    const double _a4 = 2.5 * _a2;
    const double _a5 = _a1 + _a3;
    const double _a6 = 1.0 - _es;

    // wgs-84
//    double _a = 6378137.0;
//    double _es = 6.6943799901377997e-3;
//    double _a1 = 4.2697672707157535e+4;
//    double _a2 = 1.8230912546075455e+9;
//    double _a3 = 1.4291722289812413e+2;
//    double _a4 = 4.5577281365188637e+9;
//    double _a5 = 4.2840589930055659e+4;
//    double _a6 = 9.9330562000986220e-1;

    double _x = x;
    double _y = y;
    double _z = z;

    // При необходимости переведем в Метры:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _x *= 1000.0;
            _y *= 1000.0;
            _z *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double _lon = 0;
    double _lat = 0;
    double _h = 0;

    double _zp, _w2, _w, _z2, _r2, _r, _s2, _c2, _s, _c, _ss, _g, _rg, _rf, _u, _v, _m, _f, _p;// _rm;

    _zp = std::abs( _z );
    _w2 = ( _x * _x ) + ( _y * _y );
    _w = std::sqrt( _w2 );
    _z2 = _z * _z;
    _r2 = _w2 + _z2;
    _r = std::sqrt( _r2 );
//    if( _r < 100000.0 ) {
//        lat = 0.;
//        lon = 0.;
//        h = -1.e7;
//        return;
//    }
    _lon = std::atan2( _y, _x );
    _s2 = _z2 / _r2;
    _c2 = _w2 / _r2;
    _u = _a2 / _r;
    _v = _a3 - _a4 / _r;
    if( _c2 > 0.3 ) {
        _s = ( _zp / _r ) * ( 1.0 + _c2 * ( _a1 + _u + _s2 * _v ) / _r );
        _lat = std::asin( _s );//Lat
        _ss = _s * _s;
        _c = std::sqrt( 1.0 - _ss );
    } else {
        _c = ( _w / _r ) * ( 1.0 - _s2 * ( _a5 - _u - _c2 * _v ) / _r );
        _lat = std::acos( _c );      //Lat
        _ss = 1.0 - ( _c * _c );
        _s = std::sqrt( _ss );
    }
    _g = 1.0 - ( _es * _ss );
    _rg = _a / std::sqrt( _g );
    _rf = _a6 * _rg;
    _u = _w - _rg * _c;
    _v = _zp - _rf * _s;
    _f = _c * _u + _s * _v;
    _m = _c * _v - _s * _u;
    _p = _m / ( _rf / _g + _f );
    _lat += _p; // Lat
    _h = _f + _m * _p / 2; // Height
    if( _z < 0.0 ){
        _lat = -_lat; //Lat
    }
    lat = _lat;
    lon = _lon;
    h = _h;

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat *= Convert::RdToDgD;
            lon *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            h /= 1000.0;
            break;
        }
        default:
            assert( false );
    }
}

Geodetic ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    XYZ &point )
{
    double lat, lon, h;
    ECEFtoGEO( ellipsoid, rangeUnit, angleUnit, point.X, point.Y, point.Z, lat, lon, h );
    return Geodetic( lat, lon, h );
}
//----------------------------------------------------------------------------------------------------------------------
double XYZtoDistance( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    double res = ( ( ( x1 - x2 ) * ( x1 - x2 ) ) +
                   ( ( y1 - y2 ) * ( y1 - y2 ) ) +
                   ( ( z1 - z2 ) * ( z1 - z2 ) ) );
    assert( res >= 0 );
    res = std::sqrt( res );
    return res;
}

double XYZtoDistance( const XYZ &point1, const XYZ &point2 )
{
    return XYZtoDistance( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &dX, double &dY, double &dZ )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    //double f = el.F();

    // по умолчанию Метры-Радианы:
    double _lat1 = lat1;
    double _lon1 = lon1;
    double _h1 = h1;
    double _lat2 = lat2;
    double _lon2 = lon2;
    double _h2 = h2;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat1 *= Convert::DgToRdD;
            _lon1 *= Convert::DgToRdD;
            _lat2 *= Convert::DgToRdD;
            _lon2 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h1 *= 1000.0;
            _h2 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double s1 = std::sin( lat1 );
    double c1 = std::cos( lat1 );

    double s2 = std::sin( lat2 );
    double c2 = std::cos( lat2 );

    double p1 = c1 * std::cos( lon1 );
    double p2 = c2 * std::cos( lon2 );

    double q1 = c1 * std::sin( lon1 );
    double q2 = c2 * std::sin( lon2 );

    if( Compare::AreEqualAbs( a, b ) ) { // Сфера
        dX = a * ( p2 - p1 ) + ( h2 * p2 - h1 * p1 );
        dY = a * ( q2 - q1 ) + ( h2 * q2 - h1 * q1 );
        dZ = a * ( s2 - s1 ) + ( h2 * s2 - h1 * s1 );
    } else { // Эллипсоид
        double e2 = std::pow( ellipsoid.EccentricityFirst(), 2 ); // Квадрат 1-го эксцентриситета эллипсоида

        double w1 = 1.0 / std::sqrt( 1.0 - e2 * s1 * s1 );
        double w2 = 1.0 / std::sqrt( 1.0 - e2 * s2 * s2 );

        dX = a * ( p2 * w2 - p1 * w1 ) + ( h2 * p2 - h1 * p1 );
        dY = a * ( q2 * w2 - q1 * w1 ) + ( h2 * q2 - h1 * q1 );
        dZ = ( 1.0 - e2 ) * a * ( s2 * w2 - s1 * w1 ) + ( h2 * s2 - h1 * s1 );
    }
    // dX dY dZ сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            dX *= 0.001;
            dY *= 0.001;
            dZ *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

XYZ ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point1, const Geodetic &point2 )
{
    double x, y, z;
    ECEF_offset( ellipsoid, rangeUnit, angleUnit,
        point1.Lat, point1.Lon, point1.Height, point2.Lat, point2.Lon, point2.Height, x, y, z );
    return XYZ( x, y, z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double lat, double lon, double h, double &xEast, double &yNorth, double &zUp )
{
    // по умолчанию Метры-Радианы:
    double _lat = lat;
    double _lon = lon;
    double _h = h;
    double _x = x;
    double _y = y;
    double _z = z;
    double _xr, _yr, _zr; // Reference point

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRdD;
            _lon *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h *= 1000.0;
            _x *= 1000.0;
            _y *= 1000.0;
            _z *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    GEOtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _lat, _lon, _h, _xr, _yr, _zr ); // Получены ECEF координаты опорной точки

    double cosPhi = std::cos( _lat );
    double sinPhi = std::sin( _lat );
    double cosLambda = std::cos( _lon );
    double sinLambda = std::sin( _lon );

    double _dx = _x - _xr;
    double _dy = _y - _yr;
    double _dz = _z - _zr;

    double t = ( cosLambda * _dx ) + ( sinLambda * _dy );
    xEast = ( -sinLambda * _dx ) + ( cosLambda * _dy );
    yNorth = ( -sinPhi * t ) + ( cosPhi * _dz );
    zUp = ( cosPhi * t ) + ( sinPhi * _dz );
    // xEast yNorth zUp сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            xEast *= 0.001;
            yNorth *= 0.001;
            zUp *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

ENU ECEFtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &ecef, const Geodetic &point )
{
    double e, n, u;
    ECEFtoENU( ellipsoid, rangeUnit, angleUnit, ecef.X, ecef.Y, ecef.Z, point.Lat, point.Lon, point.Height, e, n, u );
    return ENU( e, n, u );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoENUV( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double dX, double dY, double dZ, double lat, double lon, double &xEast, double &yNorth, double &zUp )
{
    // по умолчанию Метры-Радианы:
    double _lat = lat;
    double _lon = lon;
    double _dX = dX;
    double _dY = dY;
    double _dZ = dZ;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRdD;
            _lon *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _dX *= 1000.0;
            _dY *= 1000.0;
            _dZ *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double cosPhi = std::cos( _lat );
    double sinPhi = std::sin( _lat );
    double cosLambda = std::cos( _lon );
    double sinLambda = std::sin( _lon );

    double t = ( cosLambda * _dX ) + ( sinLambda * _dY );
    xEast = ( -sinLambda * _dX ) + ( cosLambda * _dY );

    zUp    =  ( cosPhi * t ) + ( sinPhi * _dZ );
    yNorth =  ( -sinPhi * t ) + ( cosPhi * _dZ );

    // xEast yNorth zUp сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            xEast *= 0.001;
            yNorth *= 0.001;
            zUp *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

ENU ECEFtoENUV( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &shift, const Geographic &point )
{
    double e, n, u;
    ECEFtoENUV( rangeUnit, angleUnit, shift.X, shift.Y, shift.Z, point.Lat, point.Lon, e, n, u );
    return ENU( e, n, u );
}
//----------------------------------------------------------------------------------------------------------------------
void ENUtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double e, double n, double u, double lat, double lon, double h, double &x, double &y, double &z )
{
    // по умолчанию Метры-Радианы:
    double _lat = lat;
    double _lon = lon;
    double _h = h;
    double _e = e;
    double _n = n;
    double _u = u;
    double _xr, _yr, _zr; // Reference point

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRdD;
            _lon *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h *= 1000.0;
            _e *= 1000.0;
            _n *= 1000.0;
            _u *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    GEOtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian,
        _lat, _lon, _h, _xr, _yr, _zr ); // Получены ECEF координаты опорной точки

    double cosPhi = std::cos( _lat );
    double sinPhi = std::sin( _lat );
    double cosLambda = std::cos( _lon );
    double sinLambda = std::sin( _lon );

    x = -sinLambda * _e - sinPhi * cosLambda * _n + cosPhi * cosLambda * _u + _xr;
    y = cosLambda * _e - sinPhi * sinLambda * _n + cosPhi * sinLambda * _u + _yr;
    z = cosPhi * _n + sinPhi * _u + _zr;

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            x *= 0.001;
            y *= 0.001;
            z *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

XYZ ENUtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &enu, const Geodetic &point )
{
    double x, y, z;
    ENUtoECEF( ellipsoid, rangeUnit, angleUnit, enu.E, enu.N, enu.U, point.Lat, point.Lon, point.Height, x, y, z );
    return XYZ( x, y, z );
}

//----------------------------------------------------------------------------------------------------------------------
void ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double &az, double &elev, double &slantRange )
{
    // по умолчанию Метры-Радианы:
    double _xEast = xEast;
    double _yNorth = yNorth;
    double _zUp = zUp;

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            _xEast *=  1000.0;
            _yNorth *= 1000.0;
            _zUp *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

//    r = std::sqrt( ( _xEast * _xEast ) + ( _yNorth * _yNorth ) ); // dangerous
    double r = std::hypot( _xEast, _yNorth ); // C++11 style

//    slantRange = sqrt( ( r * r ) + ( _zUp * _zUp ) ); // dangerous
    slantRange = std::hypot( r, _zUp ); // C++11 style
    elev = std::atan2( _zUp, r );
    az = Convert::AngleTo360( std::atan2( _xEast, _yNorth ), Units::TAngleUnit::AU_Radian );

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            slantRange *= 0.001;
            break;
        }
        default:
            assert( false );
    }
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgD;
            elev *= Convert::RdToDgD;
        }
    }
}

AER ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const ENU &point )
{
    double a, e, r;
    ENUtoAER( rangeUnit, angleUnit, point.E, point.N, point.U, a, e, r );
    return AER( a, e, r );
}
//----------------------------------------------------------------------------------------------------------------------
void AERtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double az, double elev, double slantRange, double &xEast, double &yNorth, double &zUp )
{
    double _az = az;
    double _elev = elev;
    double _slantRange = slantRange;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _az *= Convert::DgToRdD;
            _elev *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _slantRange *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    zUp = _slantRange * std::sin( _elev );
    double _r = _slantRange * std::cos( _elev );
    xEast = _r * std::sin( _az );
    yNorth = _r * std::cos( _az );
    // xEast yNorth zUp сейчас в метрах

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            xEast *= 0.001;
            yNorth *= 0.001;
            zUp *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

ENU AERtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const AER &aer )
{
    double e, n, u;
    AERtoENU( rangeUnit, angleUnit, aer.A, aer.E, aer.R, e, n, u );
    return ENU( e, n, u );
}

//----------------------------------------------------------------------------------------------------------------------
void GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double lat0, double lon0, double h0, double &xEast, double &yNorth, double &zUp )
{
    // по умолчанию Метры-Радианы:
    double _lat = lat;
    double _lon = lon;
    double _h = h;
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat *= Convert::DgToRdD;
            _lon *= Convert::DgToRdD;
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h *= 1000.0;
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    double _x, _y, _z, _x0, _y0, _z0;
    GEOtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _lat, _lon, _h, _x, _y, _z );
    GEOtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _lat0, _lon0, _h0, _x0, _y0, _z0 );

    double _dx = _x - _x0;
    double _dy = _y - _y0;
    double _dz = _z - _z0;

    ECEFtoENUV( Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _dx, _dy, _dz, _lat0, _lon0, xEast, yNorth, zUp );

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            xEast *= 0.001;
            yNorth *= 0.001;
            zUp *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

ENU GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point, const Geodetic &anchor )
{
    double e, n, u;
    GEOtoENU( ellipsoid, rangeUnit, angleUnit,
        point.Lat, point.Lon, point.Height, anchor.Lat, anchor.Lon, anchor.Height, e, n, u );
    return ENU( e, n, u );
}
//----------------------------------------------------------------------------------------------------------------------
void ENUtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double lat0, double lon0, double h0, double &lat, double &lon, double &h )
{
    // по умолчанию Метры-Радианы:
    double _xEast = xEast;
    double _yNorth = yNorth;
    double _zUp = zUp;
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _xEast *= 1000.0;
            _yNorth *= 1000.0;
            _zUp *= 1000.0;
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    double _x, _y, _z;
    ENUtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _xEast, _yNorth, _zUp,
        _lat0, _lon0, _h0, _x, _y, _z );
    ECEFtoGEO( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _x, _y, _z, lat, lon, h );

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat *= Convert::RdToDgD;
            lon *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            h *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

Geodetic ENUtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &point, const Geodetic &anchor )
{
    double lat, lon, h;
    ENUtoGEO( ellipsoid, rangeUnit, angleUnit, point.E, point.N, point.U, anchor.Lat, anchor.Lon, anchor.Height, lat, lon, h );
    return Geodetic( lat, lon, h );
}

//----------------------------------------------------------------------------------------------------------------------
void GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &az, double &elev, double &slantRange )
{
    // по умолчанию Метры-Радианы:
    double _lat1 = lat1;
    double _lon1 = lon1;
    double _h1 = h1;
    double _lat2 = lat2;
    double _lon2 = lon2;
    double _h2 = h2;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat1 *= Convert::DgToRdD;
            _lon1 *= Convert::DgToRdD;
            _lat2 *= Convert::DgToRdD;
            _lon2 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h1 *= 1000.0;
            _h2 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double _xEast, _yNorth, _zUp;
    GEOtoENU( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _lat1, _lon1, _h1, _lat2, _lon2, _h2,
        _xEast, _yNorth, _zUp );
    ENUtoAER( Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _xEast, _yNorth, _zUp, az, elev, slantRange );

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgD;
            elev *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            slantRange *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

AER GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const Geodetic &point1, const Geodetic &point2 )
{
    double a, e, r;
    GEOtoAER( ellipsoid, rangeUnit, angleUnit,
        point1.Lat, point1.Lon, point1.Height, point2.Lat, point2.Lon, point2.Height, a, e, r );
    return AER( a, e, r );
}
//----------------------------------------------------------------------------------------------------------------------
void AERtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     double az, double elev, double slantRange, double lat0, double lon0, double h0, double &lat, double &lon, double &h )
{
    // по умолчанию Метры-Радианы:
    double _az = az;
    double _elev = elev;
    double _slantRange = slantRange;
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _az *= Convert::DgToRdD;
            _elev *= Convert::DgToRdD;
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _slantRange *= 1000.0;
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double _x, _y, _z;
    AERtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _az, _elev, _slantRange,
        _lat0, _lon0, _h0, _x, _y, _z );
    ECEFtoGEO( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _x, _y, _z, lat, lon, h );

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat *= Convert::RdToDgD;
            lon *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            h *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

Geodetic AERtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const AER &aer, const Geodetic &anchor )
{
    double lat, lon, h;
    AERtoGEO( ellipsoid, rangeUnit, angleUnit, aer.A, aer.E, aer.R, anchor.Lat, anchor.Lon, anchor.Height, lat, lon, h );
    return Geodetic( lat, lon, h );
}
//----------------------------------------------------------------------------------------------------------------------
void AERtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     double az, double elev, double slantRange, double lat0, double lon0, double h0, double &x, double &y, double &z )
{
    // по умолчанию Метры-Радианы:
    double _az = az;
    double _elev = elev;
    double _slantRange = slantRange;
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _az *= Convert::DgToRdD;
            _elev *= Convert::DgToRdD;
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _slantRange *= 1000.0;
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double _x0, _y0, _z0, _e, _n, _u, _dx, _dy, _dz;
    GEOtoECEF( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _lat0, _lon0, _h0, _x0, _y0, _z0 );
    AERtoENU( Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _az, _elev, _slantRange, _e, _n, _u );
    ENUtoUVW( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _e, _n, _u, _lat0, _lon0, _dx, _dy, _dz );
    // Origin + offset from origin equals position in ECEF
    x = _x0 + _dx;
    y = _y0 + _dy;
    z = _z0 + _dz;

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            x *= 0.001;
            y *= 0.001;
            z *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

XYZ AERtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const AER &aer, const Geodetic &anchor )
{
    double x, y, z;
    AERtoECEF( ellipsoid, rangeUnit, angleUnit, aer.A, aer.E, aer.R, anchor.Lat, anchor.Lon, anchor.Height, x, y, z );
    return XYZ( x, y, z );
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double lat0, double lon0, double h0, double &az, double &elev, double &slantRange )
{
    // по умолчанию Метры-Радианы:
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;
    double _x = x;
    double _y = y;
    double _z = z;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h0 *= 1000.0;
            _x *= 1000.0;
            _y *= 1000.0;
            _z *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    double _e, _n, _u;
    ECEFtoENU( ellipsoid, Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _x, _y, _z, _lat0, _lon0, _h0,
        _e, _n, _u );
    ENUtoAER( Units::TRangeUnit::RU_Meter, Units::TAngleUnit::AU_Radian, _e, _n, _u, az, elev, slantRange );

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgD;
            elev *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            slantRange *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

AER ECEFtoAER(const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const XYZ &ecef, const Geodetic &anchor )
{
    double a, e, r;
    ECEFtoAER( ellipsoid, rangeUnit, angleUnit, ecef.X, ecef.Y, ecef.Z, anchor.Lat, anchor.Lon, anchor.Height, a, e, r );
    return AER( a, e, r );
}
//----------------------------------------------------------------------------------------------------------------------
void ENUtoUVW( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double lat0, double lon0, double &u, double &v, double &w )
{
    double _xEast = xEast;
    double _yNorth = yNorth;
    double _zUp = zUp;
    double _lat0 = lat0;
    double _lon0 = lon0;

    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _xEast *= 1000.0;
            _yNorth *= 1000.0;
            _zUp *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    double t = std::cos( _lat0 ) * _zUp - std::sin( _lat0 ) * _yNorth;
    w = std::sin( _lat0 ) * _zUp + std::cos( _lat0 ) * _yNorth;
    u = std::cos( _lon0 ) * t - std::sin( _lon0 ) * _xEast;
    v = std::sin( _lon0 ) * t + std::cos( _lon0 ) * _xEast;

    // Проверим, нужен ли перевод:
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer ):
        {
            w *= 0.001;
            u *= 0.001;
            v *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

UVW ENUtoUVW( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &enu, const Geographic &point )
{
    double u, v, w;
    ENUtoUVW( ellipsoid, rangeUnit, angleUnit, enu.E, enu.N, enu.U, point.Lat, point.Lon, u, v, w );
    return UVW( u, v, w );
}
//----------------------------------------------------------------------------------------------------------------------
double CosAngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    // Исходя из формулы косинуса угла между векторами:
    double a1 =  x1 * x2;
    double a2 =  y1 * y2;
    double a3 =  z1 * z2;
    double b1 = std::sqrt( ( x1 * x1 ) + ( y1 * y1 ) + ( z1 * z1 ) );
    double b2 = std::sqrt( ( x2 * x2 ) + ( y2 * y2 ) + ( z2 * z2 ) );

    double res = ( a1 + a2 + a3 ) / ( b1 * b2 );
    return res;
}

double CosAngleBetweenVectors( const XYZ &point1, const XYZ &point2 )
{
    return CosAngleBetweenVectors( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z );
}

//----------------------------------------------------------------------------------------------------------------------
double AngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 )
{
    return std::acos( CosAngleBetweenVectors( x1, y1, z1, x2, y2, z2 ) );
}

double AngleBetweenVectors( const XYZ &vec1, const XYZ &vec2 )
{
    return std::acos( CosAngleBetweenVectors( vec1, vec2 ) );
}
//----------------------------------------------------------------------------------------------------------------------
void VectorFromTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2, double &xV, double &yV, double &zV )
{
    xV = x2 - x1;
    yV = y2 - y1;
    zV = z2 - z1;
}

XYZ VectorFromTwoPoints( const XYZ &point1, const XYZ &point2 )
{
    XYZ result;
    VectorFromTwoPoints( point1.X, point1.Y, point1.Z, point2.X, point2.Y, point2.Z, result.X, result.Y, result.Z );
    return result;
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoECEF_3params( double xs, double ys, double zs, double dx, double dy, double dz, double &xt, double &yt, double &zt )
{
    xt = xs + dx;
    yt = ys + dy;
    zt = zs + dz;
}

void ECEFtoECEF_3params( TGeodeticDatum from, double xs, double ys, double zs, TGeodeticDatum to, double &xt, double &yt, double &zt )
{
    CShiftECEF_3 shift = GetShiftECEF_3( from, to );
    xt = xs + shift.dX();
    yt = ys + shift.dY();
    zt = zs + shift.dZ();
}

void ECEFtoECEF_3params( const TGeodeticDatum &from, XYZ ecefs, const TGeodeticDatum &to, XYZ &eceft )
{
    CShiftECEF_3 shift = GetShiftECEF_3( from, to );
    eceft.X = ecefs.X + shift.dX();
    eceft.Y = ecefs.Y + shift.dY();
    eceft.Z = ecefs.Z + shift.dZ();
}

//----------------------------------------------------------------------------------------------------------------------
CShiftECEF_3 GetShiftECEF_3( const TGeodeticDatum &from, const TGeodeticDatum &to )
{
    if( from == TGeodeticDatum::GD_SK95 && to == TGeodeticDatum::GD_PZ90 ) {
        return SK95toPZ90;
    } if( from == TGeodeticDatum::GD_PZ90 && to == TGeodeticDatum::GD_SK95 ) {
        return SK95toPZ90.Inverse();
    } else {
        assert( false );
    }
}

CShiftECEF_7 GetShiftECEF_7( const TGeodeticDatum &from, const TGeodeticDatum &to )
{
    if( from == TGeodeticDatum::GD_SK42 && to == TGeodeticDatum::GD_PZ9011 ) {
        return SK42toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_SK42 ) {
        return SK42toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_SK42 && to == TGeodeticDatum::GD_WGS84 ) {
        return SK42toWGS84;
    } else if( from == TGeodeticDatum::GD_WGS84 && to == TGeodeticDatum::GD_SK42 ) {
        return SK42toWGS84.Inverse();
    } else if( from == TGeodeticDatum::GD_SK95 && to == TGeodeticDatum::GD_PZ9011 ) {
        return SK95toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_SK95 ) {
        return SK95toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_GSK2011 && to == TGeodeticDatum::GD_PZ9011 ) {
        return GSK2011toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_GSK2011 ) {
        return GSK2011toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_PZ9002 && to == TGeodeticDatum::GD_PZ9011 ) {
        return PZ9002toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_PZ9002 ) {
        return PZ9002toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_PZ90 && to == TGeodeticDatum::GD_PZ9011 ) {
        return PZ90toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_PZ90 ) {
        return PZ90toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_WGS84 && to == TGeodeticDatum::GD_PZ9011 ) {
        return WGS84toPZ9011;
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_WGS84 ) {
        return WGS84toPZ9011.Inverse();
    } else if( from == TGeodeticDatum::GD_PZ9011 && to == TGeodeticDatum::GD_ITRF2008 ) {
        return PZ9011toITRF2008;
    } else if( from == TGeodeticDatum::GD_ITRF2008 && to == TGeodeticDatum::GD_PZ9011 ) {
        return PZ9011toITRF2008.Inverse();
    } else {
        assert( false );
    }
}
//----------------------------------------------------------------------------------------------------------------------
void ECEFtoECEF_7params( double xs, double ys, double zs, double dx, double dy, double dz,
    double rx, double ry, double rz, double s, double &xt, double &yt, double &zt )
{
    xt = ( 1.0 + s ) * ( xs - rz * ys - ry * zs ) + dx;
    yt = ( 1.0 + s ) * ( -rz * xs + ys + rx * zs ) + dy;
    zt = ( 1.0 + s ) * ( ry * xs - rx * ys + zs ) + dz;
}

void ECEFtoECEF_7params( const TGeodeticDatum &from, double xs, double ys, double zs,
    const TGeodeticDatum &to, double &xt, double &yt, double &zt )
{
    CShiftECEF_7 shift = GetShiftECEF_7( from, to );
    xt = ( 1.0 + shift.S() ) * ( xs - shift.rZ() * ys - shift.rY() * zs ) + shift.dX();
    yt = ( 1.0 + shift.S() ) * ( -shift.rZ() * xs + ys + shift.rX() * zs ) + shift.dY();
    zt = ( 1.0 + shift.S() ) * ( shift.rY() * xs - shift.rX() * ys + zs ) + shift.dZ();
};

void ECEFtoECEF_7params( const TGeodeticDatum &from, XYZ ecefs, const TGeodeticDatum &to, XYZ &eceft )
{
    CShiftECEF_7 shift = GetShiftECEF_7( from, to );
    eceft.X = ( 1.0 + shift.S() ) * ( ecefs.X - shift.rZ() * ecefs.Y - shift.rY() * ecefs.Z ) + shift.dX();
    eceft.Y = ( 1.0 + shift.S() ) * ( -shift.rZ() * ecefs.X + ecefs.Y + shift.rX() * ecefs.Y ) + shift.dY();
    eceft.Z = ( 1.0 + shift.S() ) * ( shift.rY() * ecefs.X - shift.rX() * ecefs.Y + ecefs.Z ) + shift.dZ();
}
//----------------------------------------------------------------------------------------------------------------------

void GEOtoGeoMolodenskyAbridged( const CEllipsoid &el0, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat0, double lon0, double h0, double dx, double dy, double dz,
    const CEllipsoid &el1, double &lat1, double &lon1, double &h1 )
{
    // по умолчанию Метры-Радианы:
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    double as = el0.A();
    double at = el1.A();

    double fs = 1.0 / el0.Invf();
    double ft = 1.0 / el1.Invf();

    double da = at - as;
    double df = ft - fs;

    double sinPhi = std::sin( _lat0 );
    double cosPhi = std::cos( _lat0 );
    double sinLam = std::sin( _lon0 );
    double cosLam = std::cos( _lon0 );

    double tmp = 1.0 - el0.EccentricityFirstSquared() * sinPhi * sinPhi;
    double ps = el0.A() * ( 1.0 - el0.EccentricityFirstSquared() ) / pow( tmp, 1.5 );
    double vs = el0.A() / std::sqrt( tmp );

//    double sin1sec = std::sin( 1.0 / 3600.0 * SPML::Convert::DgToRdD ); // sin of 1 sec

    // Short Molodensky formulas
//    double dlat = ( -dx * sinPhi * cosLam - dy * sinPhi * sinLam + dz * cosPhi + ( as * df + fs * da ) *
//        std::sin( 2.0 * _lat0 ) ) / ( ps * sin1sec );
//    double dlon = ( -dx * sinLam + dy * cosLam ) / ( vs * cosPhi * sin1sec );
    double dlat = ( -dx * sinPhi * cosLam - dy * sinPhi * sinLam + dz * cosPhi + ( as * df + fs * da ) *
        std::sin( 2.0 * _lat0 ) ) / ps;
    double dlon = ( -dx * sinLam + dy * cosLam ) / ( vs * cosPhi );
    double dh = dx * cosPhi * cosLam + dy * cosPhi * sinLam + dz * sinPhi + ( as * df + fs * da ) * sinPhi * sinPhi - da;

    lat1 = _lat0 + dlat;
    lon1 = _lon0 + dlon;
    h1 = _h0 + dh;

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat1 *= Convert::RdToDgD;
            lon1 *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            h1 *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}


void GEOtoGeoMolodenskyFull( const CEllipsoid &el0, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat0, double lon0, double h0, double dx, double dy, double dz, double rx, double ry, double rz, double s,
    const CEllipsoid &el1, double &lat1, double &lon1, double &h1 )
{
    // по умолчанию Метры-Радианы:
    double _lat0 = lat0;
    double _lon0 = lon0;
    double _h0 = h0;

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _lat0 *= Convert::DgToRdD;
            _lon0 *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _h0 *= 1000.0;
            break;
        }
        default:
            assert( false );
    }

    double as = el0.A();
    double at = el1.A();
    double a = ( as + at ) * 0.5;

    double fs = 1.0 / el0.Invf();
    double ft = 1.0 / el1.Invf();

    double da = at - as;
    double df = ft - fs;

    double ess = el0.EccentricityFirstSquared();
    double est = el1.EccentricityFirstSquared();
    double e2 = ( ess + est ) * 0.5;
    double de2 = est - ess;

    double sinB = std::sin( _lat0 );
    double cosB = std::cos( _lat0 );
    double cos2B = std::cos( 2.0 * _lat0 );
    double tanB = std::tan( _lat0 );
    double sinL = std::sin( _lon0 );
    double cosL = std::cos( _lon0 );

    const double ro = 180 * 60 * 60 / SPML::Consts::PI_D; // Число угловых секунд в радиане
    double tmp = 1.0 - ess * sinB * sinB;
    double M = el0.A() * ( 1.0 - ess ) / pow( tmp, 1.5 );
    double N = el0.A() / std::sqrt( tmp );

    // Full Molodensky formulas
    double dlat = ( ( N / a ) * e2 * sinB * cosB * da + ( ( N * N ) / ( a * a ) + 1.0 ) * N * sinB * cosB * de2 / 2.0 -
        ( dx * cosL + dy * sinL ) * sinB + dz * cosB ) * ro / ( M + _h0 ) -
        rx * sinL * ( 1.0 + e2 * cos2B ) + ry * cosL * ( 1.0  + e2 * cos2B ) -
        ro * s * e2 * sinB * cosB;
    double dlon = ( -dx * sinL + dy * cosL ) * ro / ( ( N + _h0 ) * cosB ) +
        tanB * ( 1.0 - e2 ) * ( rx * cosL + ry * sinL ) - rz;
    double dh = ( -a / N ) * da + N * sinB * sinB * de2 / 2.0 + ( dx * cosL + dy * sinL ) * cosB + dz * sinB -
        N * e2 * sinB * cosB * ( rx / ro * sinL - ry / ro * cosL ) + ( a * a / N  + _h0 ) * s;


    lat1 = _lat0 + dlat / 3600.0 * SPML::Convert::DgToRdD;
    lon1 = _lon0 + dlon / 3600.0 * SPML::Convert::DgToRdD;
    h1 = _h0 + dh;

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            lat1 *= Convert::RdToDgD;
            lon1 *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            h1 *= 0.001;
            break;
        }
        default:
            assert( false );
    }
}

//----------------------------------------------------------------------------------------------------------------------
void GEOtoGaussKruger( double lat, double lon, int &n, int &x, int &y )
{
    double B = lat * SPML::Convert::DgToRdD;
    double L = lon * SPML::Convert::DgToRdD;

//    n = static_cast<int>( std::ceil( ( 6.0 + lon ) / 6.0 ) );
    n = static_cast<int>( 1.0 + ( lon / 6.0 ) );
    double l = ( lon - static_cast<double>( 3 + 6 * ( n - 1 ) ) ) * SPML::Convert::DgToRdD;
    double l2 = l * l;

    double sinB = std::sin( B );
    double sin2B = std::sin( 2.0 * B );
    double sinBpow2 = sinB * sinB;
    double sinBpow4 = sinBpow2 * sinBpow2;
    double sinBpow6 = sinBpow2 * sinBpow2 * sinBpow2;

    double cosB = std::cos( B );
    double cos2B = std::cos( 2.0 * B );
    double tanB = std::tan( B );
    double sinL = std::sin( L );
    double cosL = std::cos( L );

    x = static_cast<int>(
        6367558.4968 * B - sin2B * ( 16002.8900 + 66.9607 * sinBpow2 + 0.3515 * sinBpow4 -
        l2 * ( 1594561.25 + 5336.535 * sinBpow2 + 26.790 * sinBpow4 + 0.149 * sinBpow6 +
        l2 * ( 672483.4 - 811219.9 * sinBpow2 + 5420.0 * sinBpow4 - 10.6 * sinBpow6 +
        l2 * ( 278194.0 - 830174.0 * sinBpow2 + 572434.0 * sinBpow4 - 16010.0 * sinBpow6 +
        l2 * ( 109500.0 - 574700.0 * sinBpow2 + 863700.0 * sinBpow4 - 398600.0 * sinBpow6 ))))) );

    y = static_cast<int>(
        ( 5.0 + 10.0 * n ) * 100000.0 + l * cosB * ( 6378245.0 + 21346.1415 * sinBpow2 + 107.1590 * sinBpow4 +
        0.5977 * sinBpow6 + l2 * ( 1070204.16 - 2136826.66 * sinBpow2 + 17.98 * sinBpow4 - 11.99 * sinBpow6 +
        l2 * ( 270806.0 - 1523417.0 * sinBpow2 + 1327645.0 * sinBpow4 - 21701.0 * sinBpow6 +
        l2 * ( 79690.0 - 866190.0 * sinBpow2 + 1730360.0 * sinBpow4 - 945460.0 * sinBpow6 )))) );
}

void GaussKrugerToGEO( int x, int y, double &lat, double &lon )
{
    double beta = static_cast<double>( x ) / 6367558.4968;
    double sinbeta = std::sin( beta );
    double sin2beta = std::sin( 2.0 * beta );
    double sinbetapow2 = sinbeta * sinbeta;
    double sinbetapow4 = sinbetapow2 * sinbetapow2;

    double B0 = beta + sin2beta * ( 0.00252588685 - 0.00001491860 * sinbetapow2 + 0.00000011904 * sinbetapow4 );

    double sinB0 = std::sin( B0 );
    double sin2B0 = std::sin( 2.0 * B0 );
    double sinB0pow2 = sinB0 * sinB0;
    double sinB0pow4 = sinB0pow2 * sinB0pow2;
    double sinB0pow6 = sinB0pow2 * sinB0pow2 * sinB0pow2;
    double cosB0 = std::cos( B0 );

    int n = static_cast<int>( static_cast<double>( y ) / 1000000.0 );
    double z0 = ( static_cast<double>( y ) - ( static_cast<double>( 10 * n + 5 ) * 100000.0 ) ) / ( 6378245.0 * cosB0 );
    double z02 = z0 * z0;

    double dB = -z02 * sin2B0 * ( 0.251684631 - 0.003369263 * sinB0pow2 + 0.00001127 * sinB0pow4 -
        z02 * ( 0.10500614 - 0.04559916 * sinB0pow2 + 0.00228901 * sinB0pow4 - 0.00002987 * sinB0pow6 -
        z02 * ( 0.042858 - 0.025318 * sinB0pow2 + 0.014346 * sinB0pow4 - 0.001264 * sinB0pow6 -
        z02 * ( 0.01672 - 0.00630 * sinB0pow2 + 0.01188 * sinB0pow4 - 0.00328 * sinB0pow6 ))));

    double l = z0 * ( 1.00000000000 - 0.0033467108 * sinB0pow2 - 0.0000056002 * sinB0pow4 - 0.0000000187 * sinB0pow6 -
        z02 * ( 0.16778975 + 0.16273586 * sinB0pow2 - 0.00052490 * sinB0pow4 - 0.00000846 * sinB0pow6 -
        z02 * ( 0.0420025 + 0.1487407 * sinB0pow2 + 0.0059420 * sinB0pow4 - 0.0000150 * sinB0pow6 -
        z02 * ( 0.01225 + 0.09477 * sinB0pow2 + 0.03282 * sinB0pow4 - 0.00034 * sinB0pow6 -
        z02 * ( 0.0038 + 0.0524 * sinB0pow2 + 0.0482 * sinB0pow4 - 0.0032 * sinB0pow6 )))));

    lat = ( B0 + dB ) * SPML::Convert::RdToDgD;
//    lon = ( 6.0 * ( static_cast<double>( n ) - 0.5 ) ) + l * SPML::Convert::RdToDgD;
    lon = static_cast<double>( 6 * n - 3 ) + l * SPML::Convert::RdToDgD;
}

}
}
/// \}
