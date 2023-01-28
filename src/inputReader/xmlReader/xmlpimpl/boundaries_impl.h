//
// Created by dominik on 03.12.22.
//

#pragma once

#include "inputReader/xmlReader/molsim-pskel.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"
#include <array>

namespace XMLReader {
    /**
     * @brief boundaries_pimpl implements the parser which reads the boundary conditions
     */
    class boundaries_pimpl : public boundaries_pskel {
    public:
        /**
         * @brief Function that reads the top condition
         */
        void top_boundary(const ::std::string &) override;

        /**
         * @brief Function that reads the bottom condition
         */
        void bottom_boundary(const ::std::string &) override;

        /**
         * @brief Function that reads the left condition
         */
        void left_boundary(const ::std::string &) override;

        /**
         * @brief Function that reads the right condition
         */
        void right_boundary(const ::std::string &) override;

        void front_boundary(const ::std::string &) override;

        void back_boundary(const ::std::string &) override;

        /**
         * @brief function called after reading the boundary element
         */
        void post_boundaries() override;

        /**
         * @brief function to initialize the parser
         *
         * @param cells_arg shared_ptr to LinkedCellContainer to which the conditions apply
         */
        void init(std::shared_ptr<LinkedCellStrategy> &cells_arg);

    private:
        /**
         * @brief indicator array for horizontal boundaries
         */
        std::array<double, 3> hor{0, 1, 0};

        /**
         * @brief indicator array for vertical boundaries
         */
        std::array<double, 3> vert{1, 0, 0};

        /**
         * @brief indicator array for boundaries in third dimension
         */
        std::array<double, 3> dim3{0, 0, 1};

        int first = 0;
        int second = 0;

        int third = 0;

        /**
         * @brief cells to which the boundary condition applies
         */
        std::shared_ptr<LinkedCellStrategy> cells;
    };
}
