#include "Writers.hpp"
#include <highfive/H5Easy.hpp>
#include <filesystem>

using H5Easy::File;
using H5Easy::dump;
namespace fs = std::filesystem;

void DikeDataWriter::saveData(DikeData* data, const std::string& filepath){
    File file(filepath, File::Overwrite);
    dump(file, "mesh/x", data->mesh->getx());
    dump(file, "mesh/xl", data->mesh->getxl());
    dump(file, "mesh/xr", data->mesh->getxr());

    dump(file, "ny", data->ny);
    dump(file, "yc", data->yc);
    dump(file, "yb", data->yb);
    dump(file, "width", data->width);
    dump(file, "density", data->density);
    dump(file, "pressure", data->pressure);
    dump(file, "overpressure", data->overpressure);
    dump(file, "viscosity", data->viscosity);
    dump(file, "temperature", data->temperature);
    dump(file, "time", data->time);
}