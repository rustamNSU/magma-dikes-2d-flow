#pragma once
#include <string>
#include "DikeData.hpp"


class DikeDataWriter{
    public:
        void saveData(DikeData* data, const std::string& filepath);
};