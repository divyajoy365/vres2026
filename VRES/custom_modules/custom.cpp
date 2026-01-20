/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void create_cell_types( void )
{
    if( parameters.ints.find_index("random_seed") != -1 )
    {
        SeedRandom( parameters.ints("random_seed") );
    }

    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity       = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype      = NULL;
    cell_defaults.functions.custom_cell_rule      = custom_function;

    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane          = NULL;

    cell_defaults.functions.contact_function = NULL;

    if( parameters.doubles.find_index("attachment_elastic_constant") != -1 )
    {
        cell_defaults.phenotype.mechanics.attachment_elastic_constant =
            parameters.doubles("attachment_elastic_constant");
    }

 
    initialize_cell_definitions_from_pugixml();

    // T cell
    Cell_Definition* pT = find_cell_definition( "T cell" );

    if( pT == nullptr )
    {
        std::cout << "Warning: no cell definition named 'T cell' was found.\n";
    }
    else
    {
        pT->functions.update_velocity = tcell_update_velocity;
        pT->functions.contact_function = NULL;
        pT->functions.cell_division_function  = tcell_division_function;
        pT->functions.update_phenotype = cell_proliferation_based_on_IL2;

    }

    // rod cell
    Cell_Definition* pRod = find_cell_definition( "rod" );
    if( pRod == nullptr )
    {
        std::cout << "no rod definition";
    }
    else
    {
        // ensure rods cannot move on their own accord 
        pRod->phenotype.motility.is_motile = false; 

        // set max attachments high to allow for all rod attachments including bracing attachments 
        int idx = pRod->custom_data.find_variable_index( "max_attachments" );
        if( idx < 0 )
        {
            pRod->custom_data.add_variable( "max_attachments", "dimensionless", 300);
        }
        else
        {
            pRod->custom_data[ "max_attachments" ] = 300;
        }

        // contact function enables elastic attachments between rod cells
        pRod->functions.contact_function = contact_function;

        // initialise rod cell secretion of IL-2
        pRod->functions.update_phenotype = secretion_rate_rod_cells;

    }

    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    setup_cell_rules();
    display_cell_definitions( std::cout );
    return;
}

void setup_microenvironment( void )
{
    initialize_microenvironment();
    return;
}

void setup_tissue( void )
{

    // set domain
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];

    // Working in 2D for now
    Zmin = 0.0;
    Zmax = 0.0;

    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;

    load_cells_from_pugixml();

    // ***  T cell ***

    // T cell placement, uniform random distribution across domain
    Cell_Definition* pT = find_cell_definition( "T cell" );

    if( pT != nullptr )
    {
        int number_of_t_cells = 0;
        // int number_of_t_cells = compute_number_of_tcells_from_density();
        if( parameters.ints.find_index("number_of_cells") != -1 )
        {
            number_of_t_cells = parameters.ints("number_of_cells");
        }

        for( int n = 0; n < number_of_t_cells; n++ )
        {
            // asssign random positions 
            Cell* pC = create_cell( *pT );
            double x = Xmin + UniformRandom() * Xrange;
            double y = Ymin + UniformRandom() * Yrange;
            double z = 0.0;

            pC->assign_position( x, y, z );
        }
    }
    else
    {
        std::cout << "setup_tissue error: no t cell.\n";
    }

    // ***  Rod cell ***

    // create rod cells
    Cell_Definition* pRodDef = find_cell_definition( "rod" );
    if( pRodDef == nullptr )
    {
        std::cout << "setup_tissue error: no rod cell.\n";
        return;
    }

    int number_of_rods = compute_number_of_rods_from_concentration();
    // if( parameters.ints.find_index("number_of_rods") != -1 )
    // {
    //     number_of_rods = parameters.ints("number_of_rods");
    // }


    //TODO: Make this a parameter in the XML file for easy customisation 
    int rod_cells_per_row = 27; 

    // Place rods in two staggered rows

    double rad = pRodDef->phenotype.geometry.radius;
    // dx matches elastic confluent rest length (from built in func). Ensures that
    // rod cells aren't just barely touching but introduces slight compression
    // which allows for better attachment
    //TODO: don't hard code, make parameter 
    double dx  = 2.0 * rad * 0.95238095238;    
    // vertical spacing √3 / 2  * 2 (use pythag to calculate)
    double row_gap = 0.86602540378 * dx;    

    for (int rod_idx = 0; rod_idx < number_of_rods; rod_idx++)
    {
        // assign random position for rod
        double x0 = Xmin + UniformRandom() * Xrange;
        double y0 = Ymin + UniformRandom() * Yrange;
        // working in 2D
        double z0 = 0.0; 

        // incase it's set to 3D
        if (default_microenvironment_options.simulate_2D == false)
        { z0 = Zmin + UniformRandom() * Zrange; }

        // Random rod orientation in the XY plane
        // (ux,uy) = unit vector along the rod axis
        // (vx,vy) = unit vector perpendicular to the rod axis       
        double theta = 2.0 * M_PI * UniformRandom();
        double ux = cos(theta), uy = sin(theta);
        double vx = -uy, vy = ux; 

        // each rod consist of two rows of rod cells
        std::vector<Cell*> top(rod_cells_per_row, nullptr);
        std::vector<Cell*> bot(rod_cells_per_row, nullptr);

        for (int i = 0; i < rod_cells_per_row; i++)
        {
            double offset = (i - 0.5 * (rod_cells_per_row - 1)) * dx;

            // top row
            {
                Cell* c = create_cell(*pRodDef);

                // Position = center + (offset along rod axis) + (half row_gap perpendicular)
                double x = x0 + offset * ux + 0.5 * row_gap * vx;
                double y = y0 + offset * uy + 0.5 * row_gap * vy;
                c->assign_position(x, y, z0);

                // keep track of which rod cells belong to which rods
                c->custom_data["rod_id"] = rod_idx;
                c->custom_data["row"] = 0;
                c->custom_data["col"] = i;

                top[i] = c;
            }

            // bottom row (staggered by dx/2 along axis)
            {
                Cell* c = create_cell(*pRodDef);
                double x = x0 + (offset + 0.5 * dx) * ux - 0.5 * row_gap * vx;
                double y = y0 + (offset + 0.5 * dx) * uy - 0.5 * row_gap * vy;
                c->assign_position(x, y, z0);

                c->custom_data["rod_id"] = rod_idx;
                c->custom_data["row"] = 1;
                c->custom_data["col"] = i;

                bot[i] = c;
            }
        }

        // attach along rows
        for (int i = 0; i < rod_cells_per_row - 1; i++)
        {
            attach_cells(top[i], top[i + 1]);
            attach_cells(bot[i], bot[i + 1]);
        }

        // attach diagonally to help prevent bending 
        for (int i = 0; i < rod_cells_per_row; i++)
        {
            attach_cells(top[i], bot[i]);
            if (i > 0) attach_cells(top[i], bot[i - 1]);
        }

        // Extra stiffening attachments within each row:
        // Connect each cell to the next K cells ahead in the same row.
        // This makes the rod act more like a rigid body than a floppy chain.
        // Limit K = 9. when K > 9, there is weird rod behaviour. 

        // With higher rod density, when K = 9, rods seem to break apart somewhat 
        // Reduce K 
        int K = 5;
        for (int i = 0; i < rod_cells_per_row; i++)
        {
            for (int j = i + 1; j <= std::min(rod_cells_per_row -1, i + K); j++)
            {
                attach_cells(top[i], top[j]);
                attach_cells(bot[i], bot[j]);
            }
        }

    }

    return;
}


void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
    return;
}


static inline bool get_row_col(Cell* c, int& row, int& col)
{
    int k_row = c->custom_data.find_variable_index("row");
    int k_col = c->custom_data.find_variable_index("col");

    // check that these vars exist
    if(k_row < 0 || k_col < 0) return false;

    // extract row and col value for rod cell
    row = (int) std::round(c->custom_data[k_row]);
    col = (int) std::round(c->custom_data[k_col]);
    return true;
}


static void rod_elastic_contact_with_variable_rest_length(
    Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
    // Adapted from PhysiCell's elastic_contact_with_confluent_length.
    // To keep core code untouched and enable customisation 

    // check positions are 3d 
    if( pC1->position.size() != 3 || pC2->position.size() != 3 ) return;

    // displacement vector from PC1 to PC2
    std::vector<double> displacement = pC2->position;
    displacement -= pC1->position;

    int ii = find_cell_definition_index( pC1->type );
    int jj = find_cell_definition_index( pC2->type );


    // compute adhesion strengts and attachment elastic constants for each cell
    double adhesion_ii = pC1->phenotype.mechanics.attachment_elastic_constant *
                         pC1->phenotype.mechanics.cell_adhesion_affinities[jj];
    double adhesion_jj = pC2->phenotype.mechanics.attachment_elastic_constant *
                         pC2->phenotype.mechanics.cell_adhesion_affinities[ii];

    // combine 
    double effective_k = std::sqrt( adhesion_ii * adhesion_jj );

    // default rest length as in the core code 
    double rest_length = ( p1.geometry.radius + p2.geometry.radius ) * 0.9523809523809523;


    static int rod_type = get_cell_definition("rod").type;
    if(pC1->type == rod_type && pC2->type == rod_type)
    {
        int row1, col1, row2, col2;
        if(get_row_col(pC1, row1, col1) && get_row_col(pC2, row2, col2))
        {
            // Column distance along the rod axis
            int dcol = std::abs(col1 - col2);

            // dx = spacing between adjacent rod cells
            double dx = 2.0 * p1.geometry.radius * 0.95238095238;
            if (dcol == 0)
            {
                // Cells in the same column (vertical/diagonal neighbors)
                rest_length = dx;          
            }
            else
            {

                // cells separated by dcol columns
                // rest length scales linearly with column distance
                rest_length = dcol * dx;   
            }

        }
    }

    double strength = ( norm(displacement) - rest_length ) * effective_k;
    normalize( &displacement );
    axpy( &(pC1->velocity), strength, displacement );
}



void contact_function( Cell* pMe, Phenotype& phenoMe,
                       Cell* pOther, Phenotype& phenoOther, double dt )

{
    // assign contact function for rod cells
    rod_elastic_contact_with_variable_rest_length(pMe, phenoMe, pOther, phenoOther, dt);
    return;
}



void tcell_update_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    // runs once every dt_mechanics 
    standard_update_cell_velocity( pCell, phenotype, dt );

    static int rod_type = get_cell_definition("rod").type;


    int k_touch = pCell->custom_data.find_variable_index("touching_rod");
    int k_time  = pCell->custom_data.find_variable_index("attached_time");
    int k_state = pCell->custom_data.find_variable_index("state");

    // check that these variables exist 
    if( k_touch < 0 || k_time < 0 || k_state < 0 ) return;

    // reset touching tod flag 
    pCell->custom_data[k_touch] = 0.0;
    // get all nearby cells 
    auto nearby = pCell->nearby_interacting_cells();

    for( Cell* other : nearby )
    {
        // check if it's a rod cell
        if( other->type != rod_type ) continue;

        // compute distance between cells
        double dx = pCell->position[0] - other->position[0];
        double dy = pCell->position[1] - other->position[1];
        double dz = pCell->position[2] - other->position[2];
        double d  = std::sqrt(dx*dx + dy*dy + dz*dz);

        double contact_dist = pCell->phenotype.geometry.radius + other->phenotype.geometry.radius;

        // if cells are physically touching
        if( d <= contact_dist )
        {
            // make that the t cell is toucing the rod cell in the timestep 
            pCell->custom_data[k_touch] = 1.0;
            break;
        }
    }

    // If the T cell is yet to activated
    if( pCell->custom_data[k_state] < 0.5 )
    {
        // And it is touching a rod cell 
        if( pCell->custom_data[k_touch] > 0.5 )


        {   // accumulate total contact time rod
            pCell->custom_data[k_time] += dt;

            // min required contact time required for activation 
            const double THRESH_MIN = 30;

            // if contact time exceeds threshold, activate t cell 
            if( pCell->custom_data[k_time] >= THRESH_MIN )
            {
                // set t cell to activate 
                pCell->custom_data[k_state] = 1.0;

                #pragma omp critical
                std::cout
                << "[Tcell ACTIVATED]"
                << " id =" << pCell->ID
                << " total_time =" << pCell->custom_data[k_time]
                << " t =" << PhysiCell_globals.current_time
                << std::endl;

            }
        }
    }
}


void tcell_division_function( Cell* pParent, Cell* pDaughter )
{
    // called everytime there is a cell division, includes both parent and daughter cell 
    if( pParent == nullptr || pDaughter == nullptr ) return;

    static int t_type = get_cell_definition("T cell").type;
    if( pParent->type != t_type ) return;

    int k_state = pParent->custom_data.find_variable_index("state");
    int k_time  = pParent->custom_data.find_variable_index("attached_time");
    int k_touch = pParent->custom_data.find_variable_index("touching_rod");

    if( k_state < 0 ) return;

    // If parent is active, reset only the daughter to naive
    if( pParent->custom_data[k_state] > 0.5 )
    {
        pDaughter->custom_data[k_state] = 0.0;

        if( k_time  >= 0 ) pDaughter->custom_data[k_time]  = 0.0;
        if( k_touch >= 0 ) pDaughter->custom_data[k_touch] = 0.0;

        #pragma omp critical
        std::cerr << "[DIV] parent active -> daughter naive | parent id="
                  << pParent->ID << " daughter id=" << pDaughter->ID
                  << " t=" << PhysiCell_globals.current_time << "\n";
    }
}



std::vector<std::string> coloring_function( Cell* pCell )
{
    // output = { cytoplasm, nucleus, outline, text }
    std::vector<std::string> output(4, "grey");

    // dead cells 
    if( pCell->phenotype.death.dead )
    {
        output[0] = "red";
        output[2] = "darkred";
        return output;
    }

    static int t_type   = get_cell_definition("T cell").type;
    static int rod_type = get_cell_definition("rod").type;

    // t cell 
    if( pCell->type == t_type )
    {
        int k_state = pCell->custom_data.find_variable_index("state");
        bool active = (k_state >= 0 && pCell->custom_data[k_state] > 0.5);

        if( active )
        {
            // active T cell
            output[0] = "#FF8C00";   // orange
            output[2] = "#FFFFFF";   // white outline
 
        }
        else
        {
            // naive T cell
            output[0] = "blue";
            output[2] = "darkblue";
        }
    }

    // rod cells
    else if( pCell->type == rod_type )
    {
        output[0] = "green";
        output[2] = "darkgreen";
    }

    // all other cells 
    else
    {
        output[0] = "grey";
        output[2] = "darkgrey";
    }

    return output;
}



void cell_proliferation_based_on_IL2( Cell* pCell , Phenotype& phenotype, double dt )
{
    static int t_type   = get_cell_definition("T cell").type;

    if( pCell->type == t_type )
    {
        int k_state = pCell->custom_data.find_variable_index("state");

        bool active = (k_state >= 0 && pCell->custom_data[k_state] > 0.5);
        if( active )
        {

            // extract the index for the cell proliferation cycle for the live model
            static int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
            static int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
            
            // extract the index for the IL-2 substrate
            static int IL2_index = microenvironment.find_density_index("IL-2");
            
            // extract the density of the IL-2 at the cell's position
            double IL2 = pCell->nearest_density_vector()[IL2_index];
            
            // load in the parameters of rPmax and IL initialised in the config file 
            double rPmax = parameters.doubles("rPmax");
            double IP = parameters.doubles("IP");
            
            // update the proliferation rate to be dependent on IL-2
            phenotype.cycle.data.transition_rate( cycle_start_index, cycle_end_index ) = rPmax*IL2/(IP+IL2);
        }
    }



	return;	
}



void secretion_rate_rod_cells( Cell* pCell , Phenotype& phenotype, double dt )
{
	// extract the secretion rate for rod cell 
	static int IL2_index = microenvironment.find_density_index("IL-2");
	
	// load in the parameters of rPmax and IL initialised in the config file 
	double v = parameters.doubles("v");
	double q = parameters.doubles("q");
	double rho = parameters.doubles("rho");
	
	// extract the current time iteration
	double time = PhysiCell::PhysiCell_globals.current_time;
	
	// update the secretion rate
	//std::cout << "Secretion before: "<<phenotype.secretion.secretion_rates[IL2_index] << std::endl;
	
	phenotype.secretion.secretion_rates[IL2_index] = v*q*rho*exp(-v*time);

	//std::cout << "Secretion after: "<< phenotype.secretion.secretion_rates[IL2_index] << std::endl;
	
	return;
}


int compute_number_of_rods_from_concentration()
{
    // Inputs from XML
    double conc_ug_per_mL = parameters.doubles("rod_concentration_ug_per_mL");
    int cells_per_rod     = parameters.ints("rod_cells_per_rod");
    double rho_ug_per_um3 = parameters.doubles("rho"); 

    if (conc_ug_per_mL <= 0.0 || cells_per_rod <= 0 || rho_ug_per_um3 <= 0.0)
        return 0;


    double x_min = microenvironment.mesh.bounding_box[0];
    double y_min = microenvironment.mesh.bounding_box[1];
    double z_min = microenvironment.mesh.bounding_box[2];
    double x_max = microenvironment.mesh.bounding_box[3];
    double y_max = microenvironment.mesh.bounding_box[4];
    double z_max = microenvironment.mesh.bounding_box[5];

    double dx = x_max - x_min;
    double dy = y_max - y_min;
    double dz = z_max - z_min;

    double domain_volume_um3 = dx * dy * dz;

    // Convert um^3 -> mL
    // 1 um^3 = 1e-12 mL
    double domain_volume_mL = domain_volume_um3 * 1e-12;

    // Total rod mass that should exist in the domain (ug)
    double total_mass_ug = conc_ug_per_mL * domain_volume_mL;

    // Mass of ONE rod
    // Get rod cell volume from the rod cell definition
    Cell_Definition* pRodDef = find_cell_definition("rod");
    double rod_cell_vol_um3 = pRodDef->phenotype.volume.total;

    double mass_one_rod_ug = rho_ug_per_um3 * (double)cells_per_rod * rod_cell_vol_um3;
    if (mass_one_rod_ug <= 0.0) return 0;

    int n_rods = (int) std::round(total_mass_ug / mass_one_rod_ug);
    if (n_rods < 0) n_rods = 0;

    return n_rods;
}



int compute_number_of_tcells_from_density()
{
    double cells_per_mL = parameters.doubles("tcell_concentration_cells_per_mL"); 
    if (cells_per_mL <= 0.0) return 0;

    std::cout << "\n[Tcell DEBUG] cells_per_mL = " << cells_per_mL << std::endl;

    double x_min = microenvironment.mesh.bounding_box[0];
    double y_min = microenvironment.mesh.bounding_box[1];
    double z_min = microenvironment.mesh.bounding_box[2];
    double x_max = microenvironment.mesh.bounding_box[3];
    double y_max = microenvironment.mesh.bounding_box[4];
    double z_max = microenvironment.mesh.bounding_box[5];

    double dx = x_max - x_min;
    double dy = y_max - y_min;
    double dz = z_max - z_min;

     std::cout << "[Tcell DEBUG] dx=" << dx
              << " dy=" << dy
              << " dz=" << dz << std::endl;

    double volume_um3 = dx * dy * dz;
    double volume_mL  = volume_um3 * 1e-12;

     std::cout << "[Tcell DEBUG] volume_um3=" << volume_um3
              << " volume_mL=" << volume_mL << std::endl;

    int N = (int) std::round(cells_per_mL * volume_mL);

    std::cout << "[Tcell DEBUG] t cell num3=" << N
              << std::endl;
    return std::max(0, N);
}







