#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "vector3d.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <vector>

#include "3DLsystem/Lsystem.h"
#include "2DLsystem/Lsystem2D.h"

using namespace std;

using Lines2D = vector<Line2D>;

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    // Make an image
    img::EasyImage image;
    // General
    string type = configuration["General"]["type"].as_string_or_die();
    if (type == "2DLSystem"){
        // Maak een parser aan
        LParser::LSystem2D LParser2D;
        // Maak een Lsystem - 2DLsystem
        Lsystem2D lsystem;
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        // 2DLSystem
        string inputFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        // Creëer if stream
        ifstream file(inputFile);
        // Steek in LParser
        file >> LParser2D;
        // Teken L systeem
        Lines2D lines = lsystem.drawLSystem(LParser2D, lineColor);
        image = lsystem.draw2DLines(lines, size, backgroundColor);
    }
    if (type == "Wireframe"){
        // Maak een lsystem3D
        Lsystem lystem3D;
        // Maak een Lsystem - 2DLsystem voor draw2DLines function
        Lsystem2D lsystem;
        // Maak een parser aan
        LParser::LSystem3D LParser3D;
        // Define lines
        Lines2D lines;
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        for (int i = 0; i < nrFigures; ++i) {
            string figureType = configuration["Figure" + to_string(i)]["type"].as_string_or_die();
            if (figureType == "3DLSystem") {
                // 3DLSystem
                string inputFile = configuration["Figure" + to_string(i)]["inputfile"].as_string_or_die();
                vector<double> lineColor = configuration["Figure" + to_string(i)]["color"].as_double_tuple_or_die();
                // Creëer if stream
                ifstream file(inputFile);
                // Steek in LParser
                file >> LParser3D;
                lines = lystem3D.drawLSystem(LParser3D, lineColor);
            } else {
                vector<Figure> figures = lystem3D.generateFigures(configuration);
                lines = lystem3D.doProjection(figures);
                break;
            }
        }
        // Draw image
        int size = configuration["General"]["size"].as_int_or_die();
        vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        image = lsystem.draw2DLines(lines, size, backgroundColor);
    }
    return image;
}

int main(int argc, char const* argv[])
{

    int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
