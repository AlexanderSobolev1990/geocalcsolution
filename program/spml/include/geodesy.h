//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.h
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \details    http://wiki.gis.com/wiki/index.php/Geodetic_system
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_GEODESY_H
#define SPML_GEODESY_H

// System includes:
#include <cassert>
#include <string>
#include <vector>

// SPML includes:
#include <compare.h>
#include <convert.h>
#include <units.h>

namespace SPML /// Специальная библиотека программных модулей (СБ ПМ)
{
namespace Geodesy /// Геодезические функции и функции перевода координат
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Земной эллипсоид
///
class CEllipsoid
{
public:
    ///
    /// \brief Имя эллипсоида
    /// \return Возвращает строку с именем эллипсоида
    ///
    std::string Name() const
    {
        return name;
    }

    ///
    /// \brief Большая полуось эллипсоида (экваториальный радиус)
    /// \return Возвращает большую полуось эллипсоида (экваториальный радиус) в [м]
    ///
    double A() const
    {
        return a;
    }

    ///
    /// \brief Малая полуось эллипсоида (полярный радиус)
    /// \return Возвращает малую полуось эллипсоида (полярный радиус) в [м]
    ///
    double B() const
    {
        return b;
    }

    ///
    /// \brief Cжатие f = ( a - b ) / a
    /// \return Возвращает сжатие
    ///
    double F() const
    {
        return f;
    }

    ///
    /// \brief Обратное сжатие Invf = a / ( a - b )
    /// \return Возвращает обратное сжатие
    ///
    double Invf() const
    {
        return invf;
    }

    ///
    /// \brief Первый эксцентриситет эллипсоида e1 = sqrt( ( a * a ) - ( b * b ) ) / a;
    /// \return Возвращает первый эксцентриситет эллипсоида
    ///
    double EccentricityFirst() const
    {
        return ( std::sqrt( ( a * a ) - ( b * b ) ) / a );
    }

    ///
    /// \brief Квадрат первого эксцентриситета эллипсоида es1 = 1 - ( ( b * b ) / ( a * a ) );
    /// \return Возвращает квадрат первого эксцентриситета эллипсоида
    ///
    double EccentricityFirstSquared() const
    {
        return ( 1.0 - ( ( b * b ) / ( a * a ) ) );
    }

    ///
    /// \brief Второй эксцентриситет эллипсоида e2 = sqrt( ( a * a ) - ( b * b ) ) / b;
    /// \return Возвращает второй эксцентриситет эллипсоида
    ///
    double EccentricitySecond() const
    {
        return ( std::sqrt( ( a * a ) - ( b * b ) ) / b );
    }

    ///
    /// \brief Квадрат второго эксцентриситета эллипсоида es2 = ( ( a * a ) / ( b * b ) ) - 1;
    /// \return Возвращает второй эксцентриситет эллипсоида
    ///
    double EccentricitySecondSquared() const
    {
        return ( ( ( a * a ) / ( b * b ) ) - 1.0 );
    }

    ///
    /// \brief Конструктор по умолчанию
    ///
    CEllipsoid();

    ///
    /// \brief Параметрический конструктор эллипсоида
    /// \param[in] ellipsoidName        - название эллипсоида
    /// \param[in] semiMajorAxis        - большая полуось (экваториальный радиус)
    /// \param[in] semiMinorAxis        - малая полуось (полярный радиус)
    /// \param[in] inverseFlattening    - обратное сжатие invf = a / ( a - b )
    /// \param[in] isInvfDef            - обратное сжатие задано (малая полуось расчитана из большой и обратного сжатия)
    ///
    CEllipsoid( std::string ellipsoidName, double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef );

private: // Доступ к параметрам эллипсоида после его создания не предполагается, поэтому private
    std::string name;   ///< Название эллипсоида
    double a;           ///< Большая полуось (экваториальный радиус), [м]
    double b;           ///< Малая полуось (полярный радиус) , [м]
    double invf;        ///< Обратное сжатие invf = a / ( a - b )
    double f;           ///< Сжатие f = ( a - b ) / a
};

namespace Ellipsoids /// Земные эллипсоиды
{
//------------------------------------------------------------------------------------------------------------------
//
//                                              Земные эллипсоиды:
//
//  1) Эллипсоид WGS84, https://epsg.io/7030-ellipsoid
//  2) Эллипсоид GRS80, https://epsg.io/7019-ellipsoid
//  3) Эллипсоид ПЗ-90, https://epsg.io/7054-ellipsoid
//  4) Эллипсоид Красовского, https://epsg.io/7024-ellipsoid
//  5) Сфера радиусом 6371000.0 [м], https://epsg.io/7035-ellipsoid
//  6) Сфера радиусом 6378000.0 [м]
//  7) Сфера радиусом большой полуоси эллипсоида Красовского 1940 (6378245.0 [м])
//

///
/// \brief Эллипсоид WGS84 (EPSG:7030)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257223563
///
static CEllipsoid WGS84()
{
    return CEllipsoid( "WGS84 (EPSG:7030)", 6378137.0, 0.0, 298.257223563, true );
}

///
/// \brief Эллипсоид GRS80 (EPSG:7019)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257222101
///
static CEllipsoid GRS80()
{
    return CEllipsoid( "GRS80 (EPSG:7019)", 6378137.0, 0.0, 298.257222101, true );
}

///
/// \brief Эллипсоид ПЗ-90 (EPSG:7054)
/// \details Главная полуось 6378136.0, обратное сжатие 298.257839303
///
static CEllipsoid PZ90()
{
    return CEllipsoid( "PZ90 (EPSG:7054)", 6378136.0, 0.0, 298.257839303, true );
}

///
/// \brief Эллипсоид Красовского 1940 (EPSG:7024)
/// \details Главная полуось 6378245.0, обратное сжатие 298.3
///
static CEllipsoid Krassowsky1940()
{
    return CEllipsoid( "Krasovsky1940 (EPSG:7024)", 6378245.0, 0.0, 298.3, true );
}


///
/// \brief Сфера радиусом 6371000.0 [м] (EPSG:7035)
/// \details Обратное сжатие - бесконечность
///
static CEllipsoid Sphere6371()
{
    return CEllipsoid( "Sphere 6371000.0 [м] (EPSG:7035)", 6371000.0, 6371000.0, 0.0, false );
}

///
/// \brief Сфера радиусом 6378000.0 [м]
/// \details Обратное сжатие - бесконечность
///

static CEllipsoid Sphere6378()
{
    return CEllipsoid( "Sphere 6378000.0 [м]", 6378000.0, 6378000.0, 0.0, false );
}

///
/// \brief Сфера радиусом большой полуоси эллипсоида Красовского 1940 (EPSG:7024)
/// \details Обратное сжатие - бесконечность
///
static CEllipsoid SphereKrassowsky1940()
{
    return CEllipsoid( "SphereRadiusKrasovsky1940 (EPSG:7024)", 6378245.0, 6378245.0, 0.0, false );
}

static CEllipsoid ADG66()
{
    return CEllipsoid( "Australian Geodetic Datum 1966/84 (AGD) (EPSG:?)", 6378160.0, 0.0, 298.25, true );
}

///
/// \brief Возвращает доступные предопределенные эллипсоидоы
/// \return Вектор предопределенных эллипсоидов
///
[[maybe_unused]]
static const __attribute__ ((unused)) std::vector<CEllipsoid> GetPredefinedEllipsoids()
{
    return std::vector<CEllipsoid>{        
        WGS84(),
        GRS80(),
        PZ90(),
        Krassowsky1940(),
        Sphere6371(),
        Sphere6378(),
        SphereKrassowsky1940(),
        ADG66()
    };
}

} // end namespace Ellipsoids

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Географические координаты (широта, долгота)
///
struct Geographic
{
    double Lat; ///< Широта
    double Lon; ///< Долгота

    ///
    /// \brief Конструктор по умолчанию
    ///
    Geographic() : Lat( 0.0 ), Lon( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param lat - широта
    /// \param lon - долгота
    ///
    Geographic( double lat, double lon ) : Lat( lat ), Lon( lon )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Геодезические координаты (широта, долгота, высота)
///
struct Geodetic : public Geographic
{
    double Height; ///< Высота

    ///
    /// \brief Конструктор по умолчанию
    ///
    Geodetic() : Geographic( 0.0, 0.0 ), Height( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param lat - широта
    /// \param lon - долгота
    /// \param h - высота
    ///
    Geodetic( double lat, double lon, double h ) : Geographic( lat, lon ), Height( h )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Радиолокационные координаты (расстояние по ортодроме, азимут, конечный азимут)
/// \details Имеют два азимута: Az - это начальный азимут в точке наблюдения, AzEnd - конечный азимут (по ортодроме) на дальности R
///
struct RAD
{
    double R;       /// Дальность
    double Az;      /// Азимут в точке наблюдения
    double AzEnd;   /// Азимут в точке объекта

    ///
    /// \brief Конструктор по умолчанию
    ///
    RAD() : R( 0.0 ), Az( 0.0 ), AzEnd( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param r - дальность
    /// \param az - азимут
    /// \param azEnd - конечный азимут
    ///
    RAD( double r, double az, double azEnd ) : R( r ), Az( az ), AzEnd( azEnd )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 3D декартовы ортогональные координаты (X, Y, Z)
///
struct XYZ
{
    double X; ///< X координата
    double Y; ///< Y координата
    double Z; ///< Z координата

    ///
    /// \brief Конструктор по умолчанию
    ///
    XYZ() : X( 0.0 ), Y( 0.0 ), Z( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param x - X координата
    /// \param y - Y координата
    /// \param z - Z координата
    ///
    XYZ( double x, double y, double z ) : X( x ), Y( y ), Z( z )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Координаты ENU (East-North-Up)
///
struct ENU
{
    double E; ///< East координата
    double N; ///< North координата
    double U; ///< Up координата

    ///
    /// \brief Конструктор по умолчанию
    ///
    ENU() : E( 0.0 ), N( 0.0 ), U( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param e - East координата
    /// \param n - North координата
    /// \param u - Up координата
    ///
    ENU( double e, double n, double u ) : E( e ), N( n ), U( u )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Координаты UVW
///
struct UVW
{
    double U; ///< U координата
    double V; ///< V координата
    double W; ///< W координата

    ///
    /// \brief Конструктор по умолчанию
    ///
    UVW() : U( 0.0 ), V( 0.0 ), W( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param u - U координата
    /// \param v - V координата
    /// \param w - W координата
    ///
    UVW( double u, double v, double w ) : U( u ), V( v ), W( w )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Локальные сферические координаты AER (Azimuth-Elevation-Range, Азимут-Угол места-Дальность)
///
struct AER
{
    double A; ///< Азимут
    double E; ///< Угол места
    double R; ///< Дальность

    ///
    /// \brief Конструктор по умолчанию
    ///
    AER() : A( 0.0 ), E( 0.0 ), R( 0.0 )
    {}

    ///
    /// \brief Параметрический конструктор
    /// \param a - азимут
    /// \param e - угол места
    /// \param r - дальность
    ///
    AER( double a, double e, double r ) : A( a ), E( e ), R( r )
    {}
};

//----------------------------------------------------------------------------------------------------------------------
//
//                                          Функции пересчета координат
//
[[maybe_unused]]
static double dummy_double; // Заглушка для списка параметров функций без перегрузки

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n\n  Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  latEnd    - широта конечной точки
/// \param[in]  lonEnd    - долгота конечной точки
/// \param[out] d         - расстояние между начальной и конечной точками по ортодроме
/// \param[out] az        - азимут из начальной точки на конечную
/// \param[out] azEnd     - азимут в конечной точке
///
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd = dummy_double );

///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n\n  Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] start     - начальная точка
/// \param[in] end       - конечная точка
/// \return Радиолокационные координаты (расстояние между начальной и конечной точками по ортодроме,
/// азимут из начальной точки на конечную, азимут в конечной точке)
///
RAD GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const Geographic &end );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n\n  Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  d         - расстояние между начальной и конечной точками по ортодроме
/// \param[in]  az        - азимут из начальной точки на конечную
/// \param[out] latEnd    - широта конечной точки
/// \param[out] lonEnd    - долгота конечной точки
/// \param[out] azEnd     - прямой азимут в конечной точке
///
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd = dummy_double );

///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n\n Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  start     - географические координаты начальной точки
/// \param[in]  rad       - радиолокационные координаты пути (расстояние между начальной и конечной точками по ортодроме,
///                         азимут из начальной точки на конечную)
/// \param[out] azEnd     - прямой азимут в конечной точке
/// \return Географические координаты конечной точки
///
Geographic RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geographic &start, const RAD &rad, double &azEnd = dummy_double );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет широты, долготы, высоты в декартовые геоцентрические координаты
/// \details    EPSG:9602,
///             При отсутствии высоты или расположении точки на поверхности эллипсоида, задать координату высоты h = 0.
///             Декартовые геоцентрические координаты (ECEF):
///             ось X - через пересечение гринвичского меридиана и экватора,
///             ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
///             ось Z - через северный полюс.
///             \n Источник - Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 191
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  lat       - широта точки
/// \param[in]  lon       - долгота точки
/// \param[in]  h         - высота точки над поверхностью эллипсоида
/// \param[out] x         - координата по оси Х
/// \param[out] y         - координата по оси Y
/// \param[out] z         - координата по оси Z
///
void GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double &x, double &y, double &z );

///
/// \brief Пересчет широты, долготы, высоты в декартовые геоцентрические координаты
/// \details    EPSG:9602,
///             При отсутствии высоты или расположении точки на поверхности эллипсоида, задать координату высоты h = 0.
///             Декартовые геоцентрические координаты (ECEF):
///             ось X - через пересечение гринвичского меридиана и экватора,
///             ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
///             ось Z - через северный полюс.
///             \n Источник - Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 191
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point     - геодезические координаты точки
/// \return ECEF координаты точки point
///
XYZ GEOtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет декартовых геоцентрических координат в широту, долготу, высоту
/// \details    Декартовые геоцентрические координаты (ECEF):
///             ось X - через пересечение гринвичского меридиана и экватора,
///             ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
///             ось Z - через северный полюс
///             \n Источник - Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic
///             Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476.
///             https://doi.org/10.1109/7.481290
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  x         - координата по оси Х
/// \param[in]  y         - координата по оси Y
/// \param[in]  z         - координата по оси Z
/// \param[out] lat       - широта точки
/// \param[out] lon       - долгота точки
/// \param[out] h         - высота точки над поверхностью эллипсоида
///
void ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double &lat, double &lon, double &h );

///
/// \brief Пересчет декартовых геоцентрических координат в широту, долготу, высоту
/// \details    Декартовые геоцентрические координаты (ECEF):
///             ось X - через пересечение гринвичского меридиана и экватора,
///             ось Y - через пересечение меридиана 90 [град] восточной долготы и экватора,
///             ось Z - через северный полюс
///             \n Источник - Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic
///             Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476.
///             https://doi.org/10.1109/7.481290
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  point     - точка в координатах ECEF
/// \return Геодезические (широта, долгота, высота) координаты точки point
///
Geodetic ECEFtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    XYZ &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление расстояния между точками в декартовых координатах
/// \details Вычисляется как \f$d = \sqrt{ (x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2 } \f$
/// \attention Единицы измерения выхода соответствуют единицам измерения входа
/// \param[in]  x1 - координата первой точки по оси Х
/// \param[in]  y1 - координата первой точки по оси Y
/// \param[in]  z1 - координата первой точки по оси Z
/// \param[in]  x2 - координата второй точки по оси Х
/// \param[in]  y2 - координата второй точки по оси Y
/// \param[in]  z2 - координата второй точки по оси Z
/// \return Расстояние между двумя точками в декартовых координатах
///
double XYZtoDistance( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Вычисление расстояния между точками в декартовых координатах
/// \details Вычисляется как \f$d = \sqrt{ (x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2 } \f$
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Расстояние между двумя точками в декартовых координатах
///
double XYZtoDistance( const XYZ &point1, const XYZ &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief ECEF смещение ( разница в декартовых ECEF координатах двух точек )
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  lat1      - широта 1 точки
/// \param[in]  lon1      - долгота 1 точки
/// \param[in]  h1        - высота 1 точки
/// \param[in]  lat2      - широта 2 точки
/// \param[in]  lon2      - долгота 2 точки
/// \param[in]  h2        - высота 2 точки
/// \param[out] dX        - смещение по оси X
/// \param[out] dY        - смещение по оси Y
/// \param[out] dZ        - смещение по оси Z
///
void ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &dX, double &dY, double &dZ );

///
/// \brief ECEF смещение ( разница в декартовых ECEF координатах двух точек )
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point1 - точка 1
/// \param[in] point2 - точка 2
/// \return ECEF смещение
///
XYZ ECEF_offset( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point1, const Geodetic &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ECEF координат точки в ENU относительно географических координат опорной точки (lat, lon)
/// \details https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  x         - ECEF координата X
/// \param[in]  y         - ECEF координата Y
/// \param[in]  z         - ECEF координата Z
/// \param[in]  lat       - широта опорной точки
/// \param[in]  lon       - долгота опорной точки
/// \param[in]  h         - высота опорной точки
/// \param[out] xEast     - ENU координата X (East)
/// \param[out] yNorth    - ENU координата Y (North)
/// \param[out] zUp       - ENU координата X (Up)
///
void ECEFtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double lat, double lon, double h, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод ECEF координат точки в ENU относительно  географических координат опорной точки point
/// \details https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] ecef - ECEF координаты
/// \param[in] point - опорная точка
/// \return Координаты ENU точки point
///
ENU ECEFtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &ecef, const Geodetic &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ECEF координат точки в ENU относительно географических координат (lat, lon)
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  dX        - смещение по оси X
/// \param[in]  dY        - смещение по оси Y
/// \param[in]  dZ        - смещение по оси Z
/// \param[in]  lat       - широта точки
/// \param[in]  lon       - долгота точки
/// \param[out] xEast     - ENU координата X (East)
/// \param[out] yNorth    - ENU координата Y (North)
/// \param[out] zUp       - ENU координата X (Up)
///
void ECEFtoENUV( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double dX, double dY, double dZ, double lat, double lon, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод ECEF координат точки в ENU относительно географических координат point
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] shift - смещение по декартовым осям
/// \param[in] point - точка
/// \return Координаты ENU точки point
///
ENU ECEFtoENUV( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const XYZ &shift, const Geographic &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ENU координат точки в ECEF относительно географических координат опорной точки (lat, lon)
/// \details https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  e         - East
/// \param[in]  n         - North
/// \param[in]  u         - Up
/// \param[in]  lat       - широта опорной точки
/// \param[in]  lon       - долгота опорной точки
/// \param[in]  h         - высота опорной точки
/// \param[out] x         - ECEF координата X
/// \param[out] y         - ECEF координата Y
/// \param[out] z         - ECEF координата X
///
void ENUtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double e, double n, double u, double lat, double lon, double h, double &x, double &y, double &z );

///
/// \brief Перевод ENU координат точки в ECEF относительно географических координат точки point
/// \details https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] enu - East, North, Up координаты точки
/// \param[in] point - точка
/// \return Координаты ECEF точки point
///
XYZ ENUtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &enu, const Geodetic &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ENU координат точки в AER координаты
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  xEast      - ENU координата X (East)
/// \param[in]  yNorth     - ENU координата Y (North)
/// \param[in]  zUp        - ENU координата X (Up)
/// \param[out] az         - азимут
/// \param[out] elev       - угол места
/// \param[out] slantRange - наклонная дальность
///
void ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double &az, double &elev, double &slantRange );

///
/// \brief Перевод ENU координат точки в AER координаты
/// \param[in] rangeUnit  - единицы измерения дальности
/// \param[in] angleUnit  - единицы измерения углов
/// \param[in] point      - точка в координатах ENU
/// \return Координаты AER точки point
///
AER ENUtoAER( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const ENU &point );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод AER координат точки в ENU координаты
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  az         - азимут (A)
/// \param[in]  elev       - угол места (E)
/// \param[in]  slantRange - наклонная дальность (R)
/// \param[out] xEast      - ENU координата X (East)
/// \param[out] yNorth     - ENU координата Y (North)
/// \param[out] zUp        - ENU координата X (Up)
///
void AERtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double az, double elev, double slantRange, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод ENU коордиат точки в AER координаты
/// \param[in] rangeUnit  - единицы измерения дальности
/// \param[in] angleUnit  - единицы измерения углов
/// \param[in] aer        - точка в координатах aer
/// \return Координаты ENU точки point
///
ENU AERtoENU( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit, const AER &aer );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод геодезических координат GEO точки point в координаты ENU относительно опорной точки
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  lat       - широта точки
/// \param[in]  lon       - долгота точки
/// \param[in]  h         - высота точки
/// \param[in]  lat0      - широта опорной точки
/// \param[in]  lon0      - долгота опорной точки
/// \param[in]  h0        - высота опорной точки
/// \param[out] xEast     - ENU координата X (East)
/// \param[out] yNorth    - ENU координата Y (North)
/// \param[out] zUp       - ENU координата X (Up)
///
void GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, double h, double lat0, double lon0, double h0, double &xEast, double &yNorth, double &zUp );

///
/// \brief Перевод геодезических координат GEO точки point в координаты ENU относительно опорной точки
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point     - геодезические координаты точки
/// \param[in] anchor    - геодезические координаты опорной точки, относительно которой переводим
/// \return Координаты ENU точки point
///
ENU GEOtoENU( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const Geodetic &point, const Geodetic &anchor );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод координат ENU в геодезические координаты GEO относительно опорной точки
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  xEast     - ENU координата X (East)
/// \param[in]  yNorth    - ENU координата Y (North)
/// \param[in]  zUp       - ENU координата X (Up)
/// \param[in]  lat0      - широта опорной точки
/// \param[in]  lon0      - долгота опорной точки
/// \param[in]  h0        - высота опорной точки
/// \param[out] lat       - широта точки
/// \param[out] lon       - долгота точки
/// \param[out] h         - высота точки

///
void ENUtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double lat0, double lon0, double h0, double &lat, double &lon, double &h );

///
/// \brief Перевод координат ENU в геодезические координаты GEO относительно опорной точки
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] point     - ENU координаты точки
/// \param[in] anchor    - геодезические координаты опорной точки
/// \return геодезические координаты точки point
///
Geodetic ENUtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &point, const Geodetic &anchor );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вычисление AER координат между двумя геодезическими точками
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  lat1       - широта 1 точки
/// \param[in]  lon1       - долгота 1 точки
/// \param[in]  h1         - высота 1 точки
/// \param[in]  lat2       - широта 2 точки
/// \param[in]  lon2       - долгота 2 точки
/// \param[in]  h2         - высота 2 точки
/// \param[out] az         - азимут из 1 точки на 2 точку
/// \param[out] elev       - угол места из 1 точки на 2 точку
/// \param[out] slantRange - наклонная дальность
///
void GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat1, double lon1, double h1, double lat2, double lon2, double h2, double &az, double &elev, double &slantRange );

///
/// \brief Вычисление AER координат между двумя геодезическими точками
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  point1     - геодезические координаты 1 точки
/// \param[in]  point2     - геодезические координаты 2 точки
/// \return Координаты AER между точками point1 и point2
///
AER GEOtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const Geodetic &point1, const Geodetic &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод AER координат в геодезические относительно опорной точки
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  az         - азимут из опорной точки на искомую точку
/// \param[in]  elev       - угол места из опорной точки на искомую точку
/// \param[in]  slantRange - наклонная дальность от опорной точки до искомой точки
/// \param[in]  lat0       - широта опорной точки
/// \param[in]  lon0       - долгота опорной точки
/// \param[in]  h0         - высота опорной точки
/// \param[out] lat        - широта искомой точки
/// \param[out] lon        - долгота искомой точки
/// \param[out] h          - высота искомой точки
///
void AERtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     double az, double elev, double slantRange, double lat0, double lon0, double h0, double &lat, double &lon, double &h );

///
/// \brief Перевод AER координат в геодезические относительно опорной точки
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  aer        - азимут, угол места, наклонная дальность от опорной точки на искомую
/// \param[in]  anchor     - геодезические координаты опорной точки
/// \return Геодезические координаты конечной точки координат AER относительно опорной точки
///
Geodetic AERtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const AER &aer, const Geodetic &anchor );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод AER координат относительно опорной точки в глобальные декартовые
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  az         - азимут из опорной точки на искомую точку
/// \param[in]  elev       - угол места из опорной точки на искомую точку
/// \param[in]  slantRange - наклонная дальность от опорной точки до искомой точки
/// \param[in]  lat0       - широта опорной точки
/// \param[in]  lon0       - долгота опорной точки
/// \param[in]  h0         - высота опорной точки
/// \param[out] x          - ECEF координата X
/// \param[out] y          - ECEF координата Y
/// \param[out] z          - ECEF координата X
///
void AERtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     double az, double elev, double slantRange, double lat0, double lon0, double h0, double &x, double &y, double &z );

///
/// \brief Перевод AER координат относительно опорной точки в глобальные декартовые
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  aer        - азимут, угол места, наклонная дальность
/// \param[in]  anchor     - геодезические координаты опорной точки
/// \return AER координаты точки в ECEF координатах относительно опорной точки
///
XYZ AERtoECEF( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const AER &aer, const Geodetic &anchor );
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод AER координат относительно опорной точки в глобальные декартовые
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  x          - ECEF координата X
/// \param[in]  y          - ECEF координата Y
/// \param[in]  z          - ECEF координата X
/// \param[in]  lat0       - широта опорной точки
/// \param[in]  lon0       - долгота опорной точки
/// \param[in]  h0         - высота опорной точки
/// \param[out] az         - азимут из опорной точки на искомую точку
/// \param[out] elev       - угол места из опорной точки на искомую точку
/// \param[out] slantRange - наклонная дальность от опорной точки до искомой точки
///
void ECEFtoAER( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double x, double y, double z, double lat0, double lon0, double h0, double &az, double &elev, double &slantRange );

///
/// \brief Перевод AER координат относительно опорной точки в глобальные декартовые
/// \param[in]  ellipsoid  - земной эллипсоид
/// \param[in]  rangeUnit  - единицы измерения дальности
/// \param[in]  angleUnit  - единицы измерения углов
/// \param[in]  ecef       - ECEF глобальные декратовы координаты
/// \param[in]  anchor     - геодезические координаты опорной точки
/// \return AER координаты точки в ECEF координатах относительно опорной точки
///
AER ECEFtoAER(const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
     const XYZ &ecef, const Geodetic &anchor );
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод ENU координат точки в UVW координаты
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  xEast     - East
/// \param[in]  yNorth    - North
/// \param[in]  zUp       - Up
/// \param[in]  lat0      - широта опорной точки
/// \param[in]  lon0      - долгота опорной точки
/// \param[out] u         - U координата
/// \param[out] v         - V координата
/// \param[out] w         - W координата
///
void ENUtoUVW( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double xEast, double yNorth, double zUp, double lat0, double lon0, double &u, double &v, double &w );

///
/// \brief Перевод ENU координат точки в UVW координаты
/// \details https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
/// \param[in] ellipsoid - земной эллипсоид
/// \param[in] rangeUnit - единицы измерения дальности
/// \param[in] angleUnit - единицы измерения углов
/// \param[in] enu - East, North, Up координаты точки
/// \param[in] point - точка
/// \return Координаты ECEF точки point
///
UVW ENUtoUVW( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    const ENU &enu, const Geographic &point );
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Косинус угла между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \return Косинус угла между векторами, [радиан]
///
double CosAngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Косинус угла между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Косинус угла между векторами, [радиан]
///
double CosAngleBetweenVectors( const XYZ &point1, const XYZ &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Угол между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \return Угол между векторами, [радиан]
///
double AngleBetweenVectors( double x1, double y1, double z1, double x2, double y2, double z2 );

///
/// \brief Угол между векторами в евклидовом пространстве
/// \details Предполагается, что оба вектора начинаются в точке (0, 0, 0)
/// \param[in] vec1 - вектор 1
/// \param[in] vec2 - вектор 2
/// \return Угол между векторами, [радиан]
///
double AngleBetweenVectors( const XYZ &vec1, const XYZ &vec2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Вектор из координат двух точек
/// \details Предполагается, что результирующий вектор начинается в точке (0, 0, 0)
/// \param[in] x1 - X координата 1 точки
/// \param[in] y1 - Y координата 1 точки
/// \param[in] z1 - Z координата 1 точки
/// \param[in] x2 - X координата 2 точки
/// \param[in] y2 - Y координата 2 точки
/// \param[in] z2 - Z координата 2 точки
/// \param[out] xV - X координата вектора с началом в точке 1 и концом в точке 2
/// \param[out] yV - Y координата вектора с началом в точке 1 и концом в точке 2
/// \param[out] zV - Z координата вектора с началом в точке 1 и концом в точке 2
///
void VectorFromTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2, double &xV, double &yV, double &zV );

///
/// \brief Вектор, полученный из координат двух точек
/// \details Предполагается, что результирующий вектор начинается в точке (0, 0, 0)
/// \param[in] point1 - 1 точка
/// \param[in] point2 - 2 точка
/// \return Вектор, полученный из координат двух точек
///
XYZ VectorFromTwoPoints( const XYZ &point1, const XYZ &point2 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 3-параметрическое преобразование декартовых координат из одной системы в другую
///
struct CShiftECEF_3
{
public:
    ///
    /// \brief Название смещение
    ///
    std::string Name() const
    {
        return name;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dX() const
    {
        return dx;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dY() const
    {
        return dy;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dZ() const
    {
        return dz;
    }

    CShiftECEF_3( std::string name_, double dx_, double dy_, double dz_ ) : dx{ dx_ }, dy{ dy_ }, dz{ dz_ }
    {}

    CShiftECEF_3 Inverse() const
    {
        return CShiftECEF_3( name, -dx, -dy, -dz );
    }

private:
    std::string name;
    double dx;
    double dy;
    double dz;
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 7-параметрическое преобразование декартовых координат из одной системы в другую
/// \details Также известно как преобразование Бурса-Вольфа (Bursa-Wolf)
///
struct CShiftECEF_7
{
public:
    ///
    /// \brief Название смещение
    ///
    std::string Name() const
    {
        return name;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dX() const
    {
        return dx;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dY() const
    {
        return dy;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double dZ() const
    {
        return dz;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double rX() const
    {
        return rx;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double rY() const
    {
        return ry;
    }

    ///
    /// \brief Смещение по оси X
    ///
    double rZ() const
    {
        return rz;
    }

    ///
    /// \brief S
    ///
    double S() const
    {
        return s;
    }

    CShiftECEF_7( std::string name_, double dx_, double dy_, double dz_, double rx_, double ry_, double rz_, double s_ ) :
        name{ name_ }, dx{ dx_ }, dy{ dy_ }, dz{ dz_ }, rx{ rx_ }, ry{ ry_ }, rz{ rz_ }, s{ s_ }
    {}

    CShiftECEF_7 Inverse() const
    {
        return CShiftECEF_7( name, -dx, -dy, -dz, -rx, -ry, -rz, -s );
    }

private:
    std::string name;
    double dx;
    double dy;
    double dz;
    double rx;
    double ry;
    double rz;
    double s;
};

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Геодезический датум
///
enum TGeodeticDatum : int
{
    GD_WGS84 = 0,
    GD_PZ90 = 1,
    GD_PZ9002 = 2,
    GD_PZ9011 = 3,
    GD_SK95 = 4,
    GD_SK42 = 5,
    GD_GSK2011 = 6,
    GD_ITRF2008 = 7
};

//----------------------------------------------------------------------------------------------------------------------
//                                             Параметры переводов
//---------------------------------
// 3-параметрическое преобразование
//---------------------------------
///
/// \brief SK-95 to PZ-90
/// \details ГОСТ Р 51794-2001
///
static const CShiftECEF_3 SK95toPZ90_mol( "SK95toPZ90", 25.90, -130.94, -81.76 );

///
/// \brief SK-42 to WGS-84
/// \details ГОСТ Р 51794-2001
///
static const CShiftECEF_3 SK42toWGS84_mol( "SK42toWGS84", 23.92, -141.27, -80.90 );//-80.91

///
/// \brief AGD66 to WGS-84
/// \details R.E. Deakin Department of Mathematical and Geospatial Sciences, RMIT University
/// GPO Box 2476V, MELBOURNE VIC 3001, AUSTRALIA
///
static const CShiftECEF_3 AGD66toWGS84_mol( "SK42toWGS84", -134.0, -48.0, 149.0 );

//---------------------------------
// 7-параметрическое преобразование
//---------------------------------

///
/// \brief SK-42 to WGS-84
/// \details EPSG:5044
///
static const CShiftECEF_7 SK42toWGS84( "SK42toWGS84", 23.57, -140.95, -79.8,
     0.0 / 3600.0 * SPML::Convert::DgToRdD,
    -0.35 / 3600.0 * SPML::Convert::DgToRdD,
    -0.79 / 3600.0 * SPML::Convert::DgToRdD,
    -0.00000022 );

///
/// \brief SK-42 to PZ-90.11
/// \details ГОСТ 32453-2017, Приложение А, подраздел А1
///
static const CShiftECEF_7 SK42toPZ9011( "SK42toPZ9011", 23.557, -140.844, -79.778,
    -0.00230 / 3600.0 * SPML::Convert::DgToRdD,
    -0.34646 / 3600.0 * SPML::Convert::DgToRdD,
    -0.79421 / 3600.0 * SPML::Convert::DgToRdD,
    -0.00000022800 );

///
/// \brief SK-95 to PZ-90.11
/// \details ГОСТ 32453-2017, Приложение А, подраздел А3
///
static const CShiftECEF_7 SK95toPZ9011( "SK95toPZ9011", 24.457, -130.784, -81.538,
    -0.00230 / 3600.0 * SPML::Convert::DgToRdD,
     0.00354 / 3600.0 * SPML::Convert::DgToRdD,
    -0.13421 / 3600.0 * SPML::Convert::DgToRdD,
    -0.00000022800 );

///
/// \brief GSK-2011 t PZ-90.11
/// \details EPSG:7705, ГОСТ 32453-2017, Приложение А, подраздел А5
///
static const CShiftECEF_7 GSK2011toPZ9011( "GSK2011toPZ9011", 0.0, 0.014, -0.008,
    -0.000562 / 3600.0 * SPML::Convert::DgToRdD,
     0.000019 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000053 / 3600.0 * SPML::Convert::DgToRdD,
    -0.0000000006 );

///
/// \brief PZ-90.02 to PZ-90.11
/// \details EPSG:7703, ГОСТ 32453-2017, Приложение Б, подраздел Б1
///
static const CShiftECEF_7 PZ9002toPZ9011( "PZ9002toPZ9011", -0.373, 0.186, -0.202,
    -0.00230 / 3600.0 * SPML::Convert::DgToRdD,
     0.00354 / 3600.0 * SPML::Convert::DgToRdD,
    -0.00421 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000000008 );

///
/// \brief PZ-90 to PZ-90.11
/// \details  EPSG:7703, ГОСТ 32453-2017, Приложение В, подраздел В1
///
static const CShiftECEF_7 PZ90toPZ9011( "PZ90toPZ9011", -1.443, 0.156, 0.222,
    -0.00230 / 3600.0 * SPML::Convert::DgToRdD,
     0.00354 / 3600.0 * SPML::Convert::DgToRdD,
    -0.134210 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000000228 );

///
/// \brief WGS-84 to PZ-90.11
/// \details ГОСТ 32453-2017, Приложение Г, подраздел Г1 + см. поправки в начале ГОСТа!
///
static const CShiftECEF_7 WGS84toPZ9011( "WGS84toPZ9011", -0.013, 0.106, 0.022,
    -0.00230 / 3600.0 * SPML::Convert::DgToRdD,
     0.00354 / 3600.0 * SPML::Convert::DgToRdD,
    -0.00421 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000000008 );

///
/// \brief PZ-90.11 to ITRF-2008
/// \details EPSG:7960, ГОСТ 32453-2017, Приложение Д, подраздел Д1
///
static const CShiftECEF_7 PZ9011toITRF2008( "PZ9011toITRF2008", -0.003, -0.001, 0.000,
     0.000019 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000042 / 3600.0 * SPML::Convert::DgToRdD,
     0.000002 / 3600.0 * SPML::Convert::DgToRdD,
    -0.000 );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Получить параметры перевода из СК 'from' в СК 'to'
/// \param from - СК, из которой переводят
/// \param to - СК, в которую переводят
/// \return Параметры перевода
///
CShiftECEF_3 GetShiftECEF_3( const TGeodeticDatum &from, const TGeodeticDatum &to );

///
/// \brief Получить параметры перевода из СК 'from' в СК 'to'
/// \param from - СК, из которой переводят
/// \param to - СК, в которую переводят
/// \return Параметры перевода
///
CShiftECEF_7 GetShiftECEF_7( const TGeodeticDatum &from, const TGeodeticDatum &to );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 3-параметрическое преобразование декартовых геоцентрических координат (простой сдвиг)
/// \details EPSG:9603
/// \param[in]  xs - X координата в исходной СК
/// \param[in]  ys - Y координата в исходной СК
/// \param[in]  zs - Z координата в исходной СК
/// \param[in]  dx - смещение по оси Х
/// \param[in]  dy - смещение по оси Y
/// \param[in]  dz - смещение по оси Z
/// \param[out] xt - X координата в конечной СК
/// \param[out] yt - Y координата в конечной СК
/// \param[out] zt - Z координата в конечной СК
///
void ECEFtoECEF_3params( double xs, double ys, double zs, double dx, double dy, double dz, double &xt, double &yt, double &zt );

///
/// \brief 3-параметрическое преобразование декартовых геоцентрических координат (простой сдвиг)
/// \details EPSG:9603
/// \param[in] from - исходный датум
/// \param[in]   xs - X координата в исходной СК
/// \param[in]   ys - Y координата в исходной СК
/// \param[in]   zs - Z координата в исходной СК
/// \param[in]   to - конечный датум
/// \param[out]  xt - X координата в конечной СК
/// \param[out]  yt - Y координата в конечной СК
/// \param[out]  zt - Z координата в конечной СК
///
void ECEFtoECEF_3params( const TGeodeticDatum &from, double xs, double ys, double zs,
    const TGeodeticDatum &to, double &xt, double &yt, double &zt );

///
/// \brief 3-параметрическое преобразование декартовых геоцентрических координат (простой сдвиг)
/// \details EPSG:9603
/// \param[in]   from - исходный датум
/// \param[in]  ecefs - исходные декартовы координаты в датуме 'from'
/// \param[in]     to - конечный датум
/// \param[out] eceft - конечные декартовы координаты в датуме 'to'
///
void ECEFtoECEF_3params( const TGeodeticDatum &from, XYZ ecefs, const TGeodeticDatum &to, XYZ &eceft );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief 7-параметрическое преобразование декартовых геоцентрических координат (Бурса-Вольфа)
/// \details EPSG:9606
/// \param[in]  xs - X координата в исходной СК
/// \param[in]  ys - Y координата в исходной СК
/// \param[in]  zs - Z координата в исходной СК
/// \param[in]  dx - смещение по оси Х
/// \param[in]  dy - смещение по оси Y
/// \param[in]  dz - смещение по оси Z
/// \param[in]  rx - поворот по оси Х
/// \param[in]  ry - поворот по оси Y
/// \param[in]  rz - поворот по оси Z
/// \param[in]   s - смещение (масштабный коэффициент)
/// \param[out] xt - X координата в конечной СК
/// \param[out] yt - Y координата в конечной СК
/// \param[out] zt - Z координата в конечной СК
///
void ECEFtoECEF_7params( double xs, double ys, double zs, double dx, double dy, double dz,
    double rx, double ry, double rz, double s, double &xt, double &yt, double &zt );

///
/// \brief 7-параметрическое преобразование декартовых геоцентрических координат (Бурса-Вольфа)
/// \details EPSG:9606
/// \param[in]  from - исходный датум
/// \param[in]  xs - X координата в исходной СК
/// \param[in]  ys - Y координата в исходной СК
/// \param[in]  zs - Z координата в исходной СК
/// \param[in]  to - конечный датум
/// \param[out] xt - X координата в конечной СК
/// \param[out] yt - Y координата в конечной СК
/// \param[out] zt - Z координата в конечной СК
///
void ECEFtoECEF_7params( const TGeodeticDatum &from, double xs, double ys, double zs,
    const TGeodeticDatum &to, double &xt, double &yt, double &zt );

///
/// \brief 7-параметрическое преобразование декартовых геоцентрических координат (Бурса-Вольфа)
/// \details EPSG:9606
/// \param[in]  from - исходный датум
/// \param[in]  ecefs - исходные декартовы координаты в датуме 'from'
/// \param[in]  to - конечный датум
/// \param[out] eceft - конечные декартовы координаты в датуме 'to'
///
void ECEFtoECEF_7params( const TGeodeticDatum &from, XYZ ecefs, const TGeodeticDatum &to, XYZ &eceft );

//----------------------------------------------------------------------------------------------------------------------
//                                        Преобразования Молоденского
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Сокращенное преобразование Молоденского для геодезических координат
/// \details EPSG:9605
///
void GEOtoGeoMolodenskyAbridged( const CEllipsoid &el0, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat0, double lon0, double h0, double dx, double dy, double dz,
    const CEllipsoid &el1, double &lat1, double &lon1, double &h1 );

///
/// \brief Полное преобразование Молоденского для геодезических координат
/// \details EPSG:9604
///
void GEOtoGeoMolodenskyFull( const CEllipsoid &el0, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat0, double lon0, double h0, double dx, double dy, double dz, double rx, double ry, double rz, double s,
    const CEllipsoid &el1, double &lat1, double &lon1, double &h1 );
//----------------------------------------------------------------------------------------------------------------------
//                       Геодезические координаты в плоские прямоугольные Гаусса-Крюгера
//----------------------------------------------------------------------------------------------------------------------

///
/// \brief Перевод геодезических координат из СК-42 (на эллипсоиде Красовского) в X-Y координаты Гаусса-Крюгера
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  lat - широта точки, градусы
/// \param[in]  lon - долгота точки, градусы
/// \param[out] n   - номер 6-градусной зоны (1..60)
/// \param[out] x   - вертикальная координата, метры
/// \param[out] y   - горизонтальная координата, метры
///
void SK42toGaussKruger( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double lat, double lon, int &n, int &x, int &y );

///
/// \brief Перевод X-Y координат Гаусса-Крюгера в геодезических координат из СК-42 (на эллипсоиде Красовского)
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  x   - вертикальная координата, метры
/// \param[in]  y   - горизонтальная координата, метры
/// \param[out] lat - широта точки, градусы
/// \param[out] lon - долгота точки, градусы
///
void GaussKrugerToSK42( const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    int x, int y, double &lat, double &lon );
//----------------------------------------------------------------------------------------------------------------------
} // end namespace SPML
} // end namespace Geodesy
#endif // SPML_GEODESY_H
/// \}
