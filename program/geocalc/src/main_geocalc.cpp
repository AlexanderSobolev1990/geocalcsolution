//----------------------------------------------------------------------------------------------------------------------
///
/// \file       main_geocalc.cpp
/// \brief      Консольный геодезический калькулятор
/// \date       21.12.22 - создан
/// \author     Соболев А.А.
/// \defgroup   geocalc Геодезический калькулятор
/// \details    Решение прямой и обратной геодезических задач на эллипсоиде, переводы координат
/// \addtogroup geocalc
/// \{
///

// System includes:
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <boost/program_options.hpp>

// SPML includes:
#include <spml.h>

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Возвращает строку, содержащую информацию о версии
/// \return Строка версии в формате DD-MM-YY-VV_COMMENTS, где DD - день, MM - месяц, YY - год, VV - версия, COMMENTS - комментарий(опционально)
///
static std::string GetVersion()
{
    return "GEOCALC_27.12.2022_v01_Develop";
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Печать в строку с задаваемым числом знаков после запятой
///
template <typename T>
std::string to_string_with_precision( const T a_value, const int n = 6 )
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Настройки программы
///
struct CCoordCalcSettings
{
    int Precision;                          ///< Число цифр после запятой при печати в консоль результата
    SPML::Units::TAngleUnit AngleUnit;    ///< Единицы измерения углов
    SPML::Units::TRangeUnit RangeUnit;    ///< Единицы измерения дальностей
    int EllipsoidNumber;                    ///< Эллипсоид на котором решаем геодезические задачи
    std::vector<double> Input;              ///< Входной массив
    std::string From;
    std::string To;

    ///
    /// \brief Конструктор по умолчанию
    ///
    CCoordCalcSettings()
    {
        Precision = 6;
        AngleUnit = SPML::Units::TAngleUnit::AU_Degree;
        RangeUnit = SPML::Units::TRangeUnit::RU_Kilometer;
        EllipsoidNumber = 0;
        Input.clear();
        From.clear();
        To.clear();
    }
};

//----------------------------------------------------------------------------------------------------------------------

int DetermineGeodeticDatum( std::string str, SPML::Geodesy::TGeodeticDatum &gd )
{
    if( str == "wgs84" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_WGS84;
    } else if( str == "pz90" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_PZ90;
    } else if( str == "pz9002" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_PZ9002;
    } else if( str == "pz9011" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_PZ9011;
    } else if( str == "sk95" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_SK95;
    } else if( str == "sk42" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_SK42;
    } else if( str == "gsk2011" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_GSK2011;
    } else if( str == "itrf2008" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_ITRF2008;
    } else if( str == "agd66" ) {
        gd = SPML::Geodesy::TGeodeticDatum::GD_AGD66;
    } else {
        std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief main - Основная функция
/// \param argc - количество аргументов командной строки
/// \param argv - аргументы командной строки
/// \return 0 - штатная работа, 1 - ошибка
///
int main( int argc, char *argv[] )
{
    CCoordCalcSettings settings; // Параметры приложения
    auto ellipsoids = SPML::Geodesy::Ellipsoids::GetPredefinedEllipsoids(); // Используемые эллипсоиды
    //------------------------------------------------------------------------------------------------------------------
    // Зададим параметры запуска приложения
    namespace po = boost::program_options;
    po::options_description desc( "~-= GEODETIC CALCULATOR =-~"
        "\n\nРешение геодезических задач и перевод координат (в двойной точности)"
        "\nSolve geodetic problems and convert coordinates (double precision)"
        "\n\nПараметры/Parameters", 220 ); // 220 - задает ширину строки вывода в терминал
    desc.add_options()
    // Справочные параметры:
    ( "help", "Показать эту справку и выйти/Show this text and exit" )
    ( "ver", "Показать версию и выйти/Show version and exit" )
    // Задающие параметры:
    ( "pr", po::value<int>( &settings.Precision )->default_value( settings.Precision ),
        "Число знаков после запятой при печати в косноль/Number of digits after dot while printing to console" )
    // Единицы входа дальности/углов
    ( "deg", "Вход в градусах (по умолчанию)/Input in degrees (default)" )
    ( "rad", "Вход в радианах/Input in radians" )
    ( "km", "Вход в километрах (по умолчанию)/Input in kilometers (default)" )
    ( "me", "Вход в метрах/Input in meters" )
    // Единицы выхода дальности/углов
    // На каком эллипсоиде считать
    ( "el", po::value<std::string>()->default_value( "wgs84" ), "Доступные эллипсоиды/Avaliable ellipsoids: wgs84, grs80, pz90, krasovsky1940, sphere6371, sphere6378" )
    ( "els", "Показать список доступных эллипсоидов и их параметры" )
    // Проверка
    ( "check", "Проверка решением обратной задачи/Check by solving inverse task" )    
    // Задачи:
    //------------------------------------------------------------------------------------------------------------------
    ( "geo2rad", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: LatStart LonStart LatEnd LonEnd" )
    ( "rad2geo", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: LatStart LonStart Range Azimuth" )
    //------------------------------------------------------------------------------------------------------------------
    ( "geo2ecef", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: Lat Lon Height" )
    ( "ecef2geo", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: X Y Z" )
    //------------------------------------------------------------------------------------------------------------------
    ( "ecefdist", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: X1 Y1 Z1 X2 Y2 Z2" )
    ( "ecefoffset", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: X1 Y1 Z1 X2 Y2 Z2" )
    //------------------------------------------------------------------------------------------------------------------
    ( "ecef2enu", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: X Y Z Lat0 Lon0" )
    ( "enu2ecef", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: E N U Lat0 Lon0" )
    //------------------------------------------------------------------------------------------------------------------
    ( "enu2aer", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: E N U" )
    ( "aer2enu", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: A E R" )
    //------------------------------------------------------------------------------------------------------------------
    ( "geo2enu", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: Lat Lon Height Lat0 Lon0 Height0" )
    ( "enu2geo", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: E N U Lat0 Lon0 Height0" )
    //------------------------------------------------------------------------------------------------------------------
    ( "geo2aer", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: Lat Lon Height Lat0 Lon0 Height0" )
    ( "aer2geo", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: A E R Lat0 Lon0 Height0" )
    //------------------------------------------------------------------------------------------------------------------
    ( "ecef2aer", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: X Y Z Lat0 Lon0 Height0" )
    ( "aer2ecef", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "args: A E R Lat0 Lon0 Height0" )
    //------------------------------------------------------------------------------------------------------------------
    ( "sk42toGK", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "SK-42 to Gauss-Kruger, args: Lat Lon" )
    ( "GKtosk42", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "Gauss-Kruger to SK-42, args: X Y" )
    //------------------------------------------------------------------------------------------------------------------
    ( "bw", po::value<std::vector<double>>( &settings.Input )->multitoken(),
        "Bursa-Wolf conversion (only meters in/out) (needs --from and --to keys), args: X Y Z" )
    ( "from", po::value<std::string>( &settings.From ),
        "see supported converions --list" )
    ( "to", po::value<std::string>( &settings.To ),
        "see supported converions --list" )
    ( "list", "show list of supported Bursa-Wolf conversions" )
    ;
    po::options_description cla; // Аргументы командной строки (сommand line arguments)
    cla.add( desc );
    po::variables_map vm;
//    po::store( po::command_line_parser( argc, argv ).options( cla ).run(), vm );
    po::store( po::parse_command_line( argc, argv, cla, po::command_line_style::unix_style ^ po::command_line_style::allow_short ), vm );
    po::notify( vm );
    //------------------------------------------------------------------------------------------------------------------
    // Обработаем аргументы запуска приложения
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "help" ) ) {
        std::cout << desc << std::endl;
        return EXIT_SUCCESS;
    }
    if( vm.count( "ver" ) ) {
        std::cout << SPML::GetVersion() << std::endl;
        std::cout << GetVersion() << std::endl;
        return EXIT_SUCCESS;
    }
    if( vm.count( "pr" ) ) {
        settings.Precision = vm["pr"].as<int>();
    }
    //------------------------------------------------------------------------------------------------------------------
    // Единицы углов
    if( vm.count( "deg" ) ) {
        settings.AngleUnit = SPML::Units::AU_Degree;
    }
    if( vm.count( "rad" ) ) {
        settings.AngleUnit = SPML::Units::AU_Radian;
    }
    //------------------------------------------------------------------------------------------------------------------
    // Единицы расстояния
    if( vm.count( "km" ) ) {
        settings.RangeUnit = SPML::Units::RU_Kilometer;
    }
    if( vm.count( "me" ) ) {
        settings.RangeUnit = SPML::Units::RU_Meter;
    }
    //------------------------------------------------------------------------------------------------------------------
    // Названия единиц расстояния/углов для вывода на печать
    std::string outrange;
    if( settings.RangeUnit == SPML::Units::RU_Kilometer ) {
        outrange = "km";
    } else if( settings.RangeUnit == SPML::Units::RU_Meter ) {
        outrange = "m";
    } else {
        assert( false );
    }
    std::string outangle;
    if( settings.AngleUnit == SPML::Units::AU_Degree ) {
        outangle = "deg";
    } else if( settings.AngleUnit == SPML::Units::AU_Radian ) {
        outangle = "rad";
    } else {
        assert( false );
    }
    //------------------------------------------------------------------------------------------------------------------
    // Эллипсоид
    if( vm.count( "el" ) ) {
        std::string elName = vm["el"].as<std::string>();
        if( elName == "wgs84" ) {
            settings.EllipsoidNumber = 0;
        } else if( elName == "grs80" ) {
            settings.EllipsoidNumber = 1;
        } else if( elName == "pz90" ) {
            settings.EllipsoidNumber = 2;
        } else if( elName == "krasovsky1940" ) {
            settings.EllipsoidNumber = 3;
        } else if( elName == "sphere6371" ) {
            settings.EllipsoidNumber = 4;
        } else if( elName == "sphere6378" ) {
            settings.EllipsoidNumber = 5;
        } else if( elName == "spherekrasovsky1940" ) {
            settings.EllipsoidNumber = 6;
        } else {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }
    }
    if( vm.count( "els" ) ) {
        std::string ellipsoidsString;
        for( int i = 0; i < ellipsoids.size(); i++ ) {
            ellipsoidsString += ( ellipsoids.at( i ) ).Name() +
                " a=" + std::to_string( ( ellipsoids.at( i ) ).A() )+
                " invf=" + std::to_string( ( ellipsoids.at( i ) ).Invf() );
            if( i != ellipsoids.size() - 1 ) {
                ellipsoidsString += "\n";
            }
        }
        std::cout << ellipsoidsString << std::endl;
        return EXIT_SUCCESS;
    }
    //------------------------------------------------------------------------------------------------------------------
    // Задачи:
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "geo2rad" ) ) {
        settings.Input = vm["geo2rad"].as<std::vector<double>>();
        if( settings.Input.size() != 4 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double r, az, azend;
        SPML::Geodesy::GEOtoRAD( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], settings.Input[3], r, az, azend );

        std::string result = "R[" + outrange + "] Az[" + outangle + "] AzEnd[" + outangle + "]:\n" +
            to_string_with_precision( r, settings.Precision ) + " " +
            to_string_with_precision( az, settings.Precision ) + " " +
            to_string_with_precision( azend, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double lat2, lon2, azend2;
            SPML::Geodesy::RADtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                settings.Input[0], settings.Input[1], r, az, lat2, lon2, azend2 );
            std::string result2 = "Lat[" + outangle + "] Lon" + outangle + "] AzEnd[" + outangle + "]:\n" +
                to_string_with_precision( lat2, settings.Precision ) + " " +
                to_string_with_precision( lon2, settings.Precision ) + " " +
                to_string_with_precision( azend2, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Lat[" + outangle + "] Lon" + outangle + "] AzEnd[" + outangle + "]:\n" +
                to_string_with_precision( settings.Input[2] - lat2, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[3] - lon2, settings.Precision ) + " " +
                to_string_with_precision( azend - azend2, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
        return EXIT_SUCCESS;
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "rad2geo" ) ) {
        settings.Input = vm["rad2geo"].as<std::vector<double>>();
        if( settings.Input.size() != 4 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double lat, lon, azend;
        SPML::Geodesy::RADtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], settings.Input[3], lat, lon, azend );

        std::string result = "Lat[" + outangle + "] Lon[" + outangle + "] AzEnd[" + outangle + "]:\n" +
            to_string_with_precision( lat, settings.Precision ) + " " +
            to_string_with_precision( lon, settings.Precision ) + " " +
            to_string_with_precision( azend, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double r, az, azend2;
            SPML::Geodesy::GEOtoRAD( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                settings.Input[0], settings.Input[1], lat, lon, r, az, azend2 );
            std::string result2 = "R[" + outrange + "] Az[" + outangle + "] AzEnd[" + outangle + "]:\n" +
                to_string_with_precision( r, settings.Precision ) + " " +
                to_string_with_precision( az, settings.Precision ) + " " +
                to_string_with_precision( azend2, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "R[" + outrange + "] Az[" + outangle + "] AzEnd[" + outangle + "]:\n" +
                to_string_with_precision( settings.Input[2] - r, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[3] - az, settings.Precision ) + " " +
                to_string_with_precision( azend - azend2, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
        return EXIT_SUCCESS;
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "geo2ecef" ) ) {
        settings.Input = vm["geo2ecef"].as<std::vector<double>>();
        if( settings.Input.size() != 3 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double x, y, z;
        SPML::Geodesy::GEOtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
             settings.Input[0], settings.Input[1], settings.Input[2], x, y, z );

        std::string result = "X[" + outangle + "] Y[" + outangle + "] Z[" + outangle + "]:\n" +
            to_string_with_precision( x, settings.Precision ) + " " +
            to_string_with_precision( y, settings.Precision ) + " " +
            to_string_with_precision( z, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double lat, lon, h;
            SPML::Geodesy::ECEFtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                x, y, z, lat, lon, h );
            std::string result2 = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( lat, settings.Precision ) + " " +
                to_string_with_precision( lon, settings.Precision ) + " " +
                to_string_with_precision( h, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - lat, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - lon, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - h, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "ecef2geo" ) ) {
        settings.Input = vm["ecef2geo"].as<std::vector<double>>();
        if( settings.Input.size() != 3 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double lat, lon, h;
        SPML::Geodesy::ECEFtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], lat, lon, h );
        std::string result = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
            to_string_with_precision( lat, settings.Precision ) + " " +
            to_string_with_precision( lon, settings.Precision ) + " " +
            to_string_with_precision( h, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double x, y, z;
            SPML::Geodesy::GEOtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                 lat, lon, h, x, y, z );
            std::string result2 = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( x, settings.Precision ) + " " +
                to_string_with_precision( y, settings.Precision ) + " " +
                to_string_with_precision( z, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - x, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - y, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - z, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "ecefdist" ) ) {
        settings.Input = vm["ecefdist"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double d = SPML::Geodesy::XYZtoDistance( settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5] );
        if( settings.RangeUnit == SPML::Units::RU_Kilometer ) {
            d *= 0.001;
        }
        std::string result = "Distance[" + outrange + "]:\n" +
            to_string_with_precision( d, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nNo check provided for this operation!" << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "ecefoffset" ) ) {
        settings.Input = vm["ecefoffset"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double dx, dy, dz;
        SPML::Geodesy::ECEF_offset( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], dx, dy, dz );
        if( settings.RangeUnit == SPML::Units::RU_Kilometer ) {
            dx *= 0.001;
            dy *= 0.001;
            dz *= 0.001;
        }
        std::string result = "dX[" + outrange + "] dY[" + outrange + "] dZ[" + outrange + "]:\n" +
            to_string_with_precision( dx, settings.Precision ) + " " +
            to_string_with_precision( dy, settings.Precision ) + " " +
            to_string_with_precision( dz, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nNo check provided for this operation!" << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------    
    if( vm.count( "ecef2enu" ) ) {
        settings.Input = vm["ecef2enu"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double e, n, u;
        SPML::Geodesy::ECEFtoENU( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], settings.Input[3], settings.Input[4], settings.Input[5],
            e, n, u );
        std::string result = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( n, settings.Precision ) + " " +
            to_string_with_precision( u, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double x, y, z;
            SPML::Geodesy::ENUtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                 e, n, u, settings.Input[3], settings.Input[4], settings.Input[5], x, y, z );
            std::string result2 = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( x, settings.Precision ) + " " +
                to_string_with_precision( y, settings.Precision ) + " " +
                to_string_with_precision( z, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - x, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - y, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - z, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "enu2ecef" ) ) {
        settings.Input = vm["enu2ecef"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double x, y, z;
        SPML::Geodesy::ENUtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], x, y, z );
        std::string result = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
            to_string_with_precision( x, settings.Precision ) + " " +
            to_string_with_precision( y, settings.Precision ) + " " +
            to_string_with_precision( z, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double e, n, u;
            SPML::Geodesy::ECEFtoENU( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                 x, y, z, settings.Input[3], settings.Input[4], settings.Input[5], e, n, u );
            std::string result2 = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( e, settings.Precision ) + " " +
                to_string_with_precision( n, settings.Precision ) + " " +
                to_string_with_precision( u, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - e, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - n, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - u, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "enu2aer" ) ) {
        settings.Input = vm["enu2aer"].as<std::vector<double>>();
        if( settings.Input.size() != 3 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double a, e, r;
        SPML::Geodesy::ENUtoAER( settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], a, e, r );
        std::string result = "Azimuth[" + outangle + "]. Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
            to_string_with_precision( a, settings.Precision ) + " " +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( r, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double e_, n_, u_;
            SPML::Geodesy::AERtoENU( settings.RangeUnit, settings.AngleUnit, a, e, r, e_, n_, u_ );
            std::string result2 = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( e_, settings.Precision ) + " " +
                to_string_with_precision( n_, settings.Precision ) + " " +
                to_string_with_precision( u_, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - e_, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - n_, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - u_, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "aer2enu" ) ) {
        settings.Input = vm["aer2enu"].as<std::vector<double>>();
        if( settings.Input.size() != 3 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double e, n, u;
        SPML::Geodesy::AERtoENU( settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], e, n, u );
        std::string result = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( n, settings.Precision ) + " " +
            to_string_with_precision( u, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double a_, e_, r_;
            SPML::Geodesy::ENUtoAER( settings.RangeUnit, settings.AngleUnit, e, n, u, a_, e_, r_ );
            std::string result2 = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
                to_string_with_precision( a_, settings.Precision ) + " " +
                to_string_with_precision( e_, settings.Precision ) + " " +
                to_string_with_precision( r_, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - a_, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - e_, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - r_, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "geo2enu" ) ) {
        settings.Input = vm["geo2enu"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double e, n, u;
        SPML::Geodesy::GEOtoENU( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], e, n, u );
        std::string result = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( n, settings.Precision ) + " " +
            to_string_with_precision( u, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double lat, lon, h;
            SPML::Geodesy::ENUtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                e, n, u, settings.Input[3], settings.Input[4], settings.Input[5], lat, lon, h );
            std::string result2 = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( lat, settings.Precision ) + " " +
                to_string_with_precision( lon, settings.Precision ) + " " +
                to_string_with_precision( h, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - lat, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - lon, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - h, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "enu2geo" ) ) {
        settings.Input = vm["enu2geo"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double lat, lon, h;
        SPML::Geodesy::ENUtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], lat, lon, h );
        std::string result = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
            to_string_with_precision( lat, settings.Precision ) + " " +
            to_string_with_precision( lon, settings.Precision ) + " " +
            to_string_with_precision( h, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double e, n, u;
            SPML::Geodesy::GEOtoENU( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                lat, lon, h, settings.Input[3], settings.Input[4], settings.Input[5], e, n, u );
            std::string result2 = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( e, settings.Precision ) + " " +
                to_string_with_precision( n, settings.Precision ) + " " +
                to_string_with_precision( u, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - e, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - n, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - u, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "geo2aer" ) ) {
        settings.Input = vm["geo2aer"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double a, e, r;
        SPML::Geodesy::GEOtoAER( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], a, e, r );
        std::string result = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
            to_string_with_precision( a, settings.Precision ) + " " +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( r, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double lat, lon, h;
            SPML::Geodesy::AERtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                a, e, r, settings.Input[3], settings.Input[4], settings.Input[5], lat, lon, h );
            std::string result2 = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( lat, settings.Precision ) + " " +
                to_string_with_precision( lon, settings.Precision ) + " " +
                to_string_with_precision( h, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - lat, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - lon, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - h, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "aer2geo" ) ) {
        settings.Input = vm["aer2geo"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double lat, lon, h;
        SPML::Geodesy::ENUtoGEO( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], lat, lon, h );
        std::string result = "Lat[" + outangle + "] Lon[" + outangle + "] Height[" + outrange + "]:\n" +
            to_string_with_precision( lat, settings.Precision ) + " " +
            to_string_with_precision( lon, settings.Precision ) + " " +
            to_string_with_precision( h, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double e, n, u;
            SPML::Geodesy::GEOtoENU( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                lat, lon, h, settings.Input[3], settings.Input[4], settings.Input[5], e, n, u );
            std::string result2 = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( e, settings.Precision ) + " " +
                to_string_with_precision( n, settings.Precision ) + " " +
                to_string_with_precision( u, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "East[" + outrange + "] North[" + outrange + "] Up[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - e, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - n, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - u, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "ecef2aer" ) ) {
        settings.Input = vm["ecef2aer"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double a, e, r;
        SPML::Geodesy::ECEFtoAER( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2], settings.Input[3], settings.Input[4], settings.Input[5],
            a, e, r );
        std::string result = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
            to_string_with_precision( a, settings.Precision ) + " " +
            to_string_with_precision( e, settings.Precision ) + " " +
            to_string_with_precision( r, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double x, y, z;
            SPML::Geodesy::AERtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                 a, e, r, settings.Input[3], settings.Input[4], settings.Input[5], x, y, z );
            std::string result2 = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( x, settings.Precision ) + " " +
                to_string_with_precision( y, settings.Precision ) + " " +
                to_string_with_precision( z, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - x, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - y, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - z, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "aer2ecef" ) ) {
        settings.Input = vm["aer2ecef"].as<std::vector<double>>();
        if( settings.Input.size() != 6 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }

        double x, y, z;
        SPML::Geodesy::AERtoECEF( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], settings.Input[2],
            settings.Input[3], settings.Input[4], settings.Input[5], x, y, z );
        std::string result = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
            to_string_with_precision( x, settings.Precision ) + " " +
            to_string_with_precision( y, settings.Precision ) + " " +
            to_string_with_precision( z, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double a, e, r;
            SPML::Geodesy::ECEFtoAER( ellipsoids.at( settings.EllipsoidNumber ), settings.RangeUnit, settings.AngleUnit,
                 x, y, z, settings.Input[3], settings.Input[4], settings.Input[5], a, e, r );
            std::string result2 = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
                to_string_with_precision( a, settings.Precision ) + " " +
                to_string_with_precision( e, settings.Precision ) + " " +
                to_string_with_precision( r, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Azimuth[" + outangle + "] Elevation[" + outangle + "] slantRange[" + outrange + "]:\n" +
                to_string_with_precision( settings.Input[0] - a, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - e, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[2] - r, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "sk42toGK" ) ) {
        settings.Input = vm["sk42toGK"].as<std::vector<double>>();
        if( settings.Input.size() != 2 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }
        int n, x, y;
        SPML::Geodesy::SK42toGaussKruger( settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], n, x, y );
        std::string result = "N X[" + outrange + "] Y[" + outrange + "]: \n" +
            std::to_string( n ) + " " +
            to_string_with_precision( x, settings.Precision ) + " " +
            to_string_with_precision( y, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            double lat, lon;
            SPML::Geodesy::GaussKrugerToSK42( settings.RangeUnit, settings.AngleUnit, x, y, lat, lon );
            std::string result2 = "Lat[" + outangle + "] Lon[" + outangle + "]:\n" +
                to_string_with_precision( lat, settings.Precision ) + " " +
                to_string_with_precision( lon, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "Lat[" + outangle + "] Lon[" + outangle + "]:\n" +
                to_string_with_precision( settings.Input[0] - lat, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - lon, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "GKtosk42" ) ) {
        settings.Input = vm["GKtosk42"].as<std::vector<double>>();
        if( settings.Input.size() != 2 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }
        double lat, lon;
        SPML::Geodesy::GaussKrugerToSK42( settings.RangeUnit, settings.AngleUnit,
            settings.Input[0], settings.Input[1], lat, lon );
        std::string result = "Lat[" + outangle + "] Lon[" + outangle + "]:\n" +
            to_string_with_precision( lat, settings.Precision ) + " " +
            to_string_with_precision( lon, settings.Precision );
        std::cout << result << std::endl;

        if( vm.count( "check" ) ) {
            std::cout << "\nCheck by solving inverse task and calc delta:" << std::endl;
            int n, x, y;
            SPML::Geodesy::SK42toGaussKruger( settings.RangeUnit, settings.AngleUnit, lat, lon, n, x, y );
            std::string result2 = "N X[" + outrange + "] Y[" + outrange + "]: \n" +
                std::to_string( n ) + " " +
                to_string_with_precision( x, settings.Precision ) + " " +
                to_string_with_precision( y, settings.Precision );
            std::cout << result2 << std::endl;
            std::cout << "\nDelta:" << std::endl;
            std::string resultDelta = "X[" + outrange + "] Y[" + outrange + "]: \n" +
                to_string_with_precision( settings.Input[0] - x, settings.Precision ) + " " +
                to_string_with_precision( settings.Input[1] - y, settings.Precision );
            std::cout << resultDelta << std::endl;
        }
    }
    SPML::Geodesy::TGeodeticDatum _from;
    SPML::Geodesy::TGeodeticDatum _to;
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "from" ) ) {
        std::string from = vm["from"].as<std::string>();
        DetermineGeodeticDatum( from, _from );
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "to" ) ) {
        std::string to = vm["to"].as<std::string>();
        DetermineGeodeticDatum( to, _to );
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "list" ) ) {
        std::cout << "Bursa-Wolf conversions:\n" <<
        "sk42 <--> wgs84 \n" <<
        "sk42 <--> pz9011 \n" <<
        "sk95 <--> pz9011 \n" <<
        "gsk2011 <--> pz9011 \n" <<
        "pz9002 <--> pz9011 \n" <<
        "pz90 <--> pz9011 \n" <<
        "wgs84 <--> pz9011 \n" <<
        "pz9011 <--> itrf2008"<< std::endl;
        return EXIT_SUCCESS;
    }
    //------------------------------------------------------------------------------------------------------------------
    if( vm.count( "bw" ) ) {
        settings.Input = vm["bw"].as<std::vector<double>>();
        if( settings.Input.size() != 3 ) {
            std::cout << "Неверный ввод, смотри --help/Wrong input, read --help" << std::endl;
            return EXIT_FAILURE;
        }
        double x, y, z;
        SPML::Geodesy::ECEFtoECEF_7params( _from, settings.Input[0], settings.Input[1], settings.Input[2], _to, x, y, z );
//        std::string result2 = "X[" + outrange + "] Y[" + outrange + "] Z[" + outrange + "]:\n" +
        std::string result = "X[m] Y[m] Z[m]:\n" +
            to_string_with_precision( x, settings.Precision ) + " " +
            to_string_with_precision( y, settings.Precision ) + " " +
            to_string_with_precision( z, settings.Precision );
        std::cout << result << std::endl;
    }
    //------------------------------------------------------------------------------------------------------------------
    return EXIT_SUCCESS;
}// end main
/// \}
