//
// Created by lukas on 01.12.22.
//
#pragma once

#include "inputReader/xmlReader/xmlpimpl/molecular_pimpl.h"
#include "inputReader/FileReader.h"
#include "inputReader/xmlReader/xmlpimpl/from_checkpoint.h"

/**
 * @brief XMLReader contains functionality for reading an xml-file
 */
namespace XMLReader {

    /**
     * @brief XmlReader combines all readers into one reader for xml-files
     */
    class XmlReader : public inputReader::FileReader {

    private:
        /**
         * @brief xml-reader for ints
         */
        xml_schema::int_pimpl int_p;

        /**
         * @brief xml-reader for strings
         */
        xml_schema::string_pimpl string_p;

        /**
         * @brief xml-reader for doubles
         */
        xml_schema::double_pimpl double_p;

        xml_schema::boolean_pimpl bool_p;

        /**
         * @brief xml-reader for simualtion
         */
        simulation_pimpl sim_p;

        /**
         * @brief xml-reader for cuboids
         */
        cuboid_pimpl cub_p;

        /**
         * @brief xml-reader for input paths to .txt files for generating cuboids
         */
        cuboid_input_pimpl cub_in_p;

        /**
         * @brief xml-reader for input paths to .txt files for generating spheres
         */
        spheres_input_pimpl sph_in_p;

        /**
         * @brief xml-reader for reading boundaries
         */
        boundaries_pimpl bou_in_p;

        /**
         * @brief xml-readr for reading the specifications for the thermostat
         */
        temperature_pimpl temp_p;

        /**
         * @brief xml-parser for reading the specifications for a sphere from xml
         */
        sphere_pimpl sph_p;

        /**
         * @brief xml-parser to read particles from a checkpoint file
         */
        from_checkpoint_pimpl fr_check;

        membrane_pimpl memb_p;

        statistics_impl stat_p;

    public:
        /**
         * @brief constructor of XmlReader
         * @param s filename of the input file
         */
        explicit XmlReader(std::string &s);

        /**
         * @brief overridden function of InputReader
         * @param sim simulation which is initialized using the input
         */
        void read(std::shared_ptr<Simulation> &sim) override;

        /**
         * @brief reads input from file and initializes the given Simulation and LinkedCellContainer
         * @param sim lvalue reference to shared_ptr pointing to simulation
         * @param lc lvalue reference to shared_ptr pointing to LinkedCellContainer
         */
        void read(std::shared_ptr<Simulation> &sim, std::shared_ptr<LinkedCellContainer> &lc);

    };

}