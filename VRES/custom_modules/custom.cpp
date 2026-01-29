// /*
// ###############################################################################
// # If you use PhysiCell in your project, please cite PhysiCell and the version #
// # number, such as below:                                                      #
// #                                                                             #
// # We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
// #                                                                             #
// # [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
// #     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
// #     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
// #     DOI: 10.1371/journal.pcbi.1005991                                       #
// #                                                                             #
// # See VERSION.txt or call get_PhysiCell_version() to get the current version  #
// #     x.y.z. Call display_citations() to get detailed information on all cite-#
// #     able software used in your PhysiCell application.                       #
// #                                                                             #
// # Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
// #     as below:                                                               #
// #                                                                             #
// # We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
// # with BioFVM [2] to solve the transport equations.                           #
// #                                                                             #
// # [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
// #     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
// #     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
// #     DOI: 10.1371/journal.pcbi.1005991                                       #
// #                                                                             #
// # [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
// #     llelized diffusive transport solver for 3-D biological simulations,     #
// #     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
// #                                                                             #
// ###############################################################################
// #                                                                             #
// # BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
// #                                                                             #
// # Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
// # All rights reserved.                                                        #
// #                                                                             #
// # Redistribution and use in source and binary forms, with or without          #
// # modification, are permitted provided that the following conditions are met: #
// #                                                                             #
// # 1. Redistributions of source code must retain the above copyright notice,   #
// # this list of conditions and the following disclaimer.                       #
// #                                                                             #
// # 2. Redistributions in binary form must reproduce the above copyright        #
// # notice, this list of conditions and the following disclaimer in the         #
// # documentation and/or other materials provided with the distribution.        #
// #                                                                             #
// # 3. Neither the name of the copyright holder nor the names of its            #
// # contributors may be used to endorse or promote products derived from this   #
// # software without specific prior written permission.                         #
// #                                                                             #
// # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
// # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
// # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
// # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
// # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
// # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
// # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
// # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
// # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
// # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
// # POSSIBILITY OF SUCH DAMAGE.                                                 #
// #                                                                             #
// ###############################################################################
// */


#include "./custom.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cstdint>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void contact_function( Cell* pMe, Phenotype& phenoMe,
                       Cell* pOther, Phenotype& phenoOther, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );
void tcell_update_velocity( Cell* pCell, Phenotype& phenotype, double dt );
void tcell_division_function( Cell* pParent, Cell* pDaughter );
std::vector<std::string> coloring_function( Cell* pCell );
void cell_proliferation_based_on_IL2( Cell* pCell , Phenotype& phenotype, double dt );
void secretion_rate_rod_cells( Cell* pCell , Phenotype& phenotype, double dt );
int compute_number_of_rods_from_concentration();
int compute_number_of_tcells_from_density();
void avoid_boundaries( Cell* pCell , Phenotype& phenotype, double dt );
void clamp_cell_to_domain(Cell* c);



// helper functions 
struct Vec3 { double x,y,z; };
static inline Vec3 operator+(Vec3 a, Vec3 b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(Vec3 a, Vec3 b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(double s, Vec3 a){ return {s*a.x,s*a.y,s*a.z}; }
static inline double dot3(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 cross3(Vec3 a, Vec3 b){
    return {a.y*b.z - a.z*b.y,
            a.z*b.x - a.x*b.z,
            a.x*b.y - a.y*b.x};
}

static inline double norm3(Vec3 a){ return std::sqrt(std::max(0.0, dot3(a,a))); }
static inline Vec3 normalize3(Vec3 a){
    double n = norm3(a);
    if(n < 1e-12) return {1,0,0};
    return (1.0/n) * a;
}


// check if two cells are attached to each other 
static inline bool are_attached(Cell* a, Cell* b)
{
    for(Cell* x : a->state.attached_cells)
        if(x == b) return true;
    return false;
}


// pre attachment rest length storage 
// This allows for customisable rest length for the elastic attachment function 

// store the rest lengths between 2 cells in the 
static std::unordered_map<std::uint64_t, double> g_rest_length;


static inline std::uint64_t pair_key(std::uint32_t a, std::uint32_t b)
{
    // use two 32-bit ids to create a 64-bit key 
    if(a > b) std::swap(a,b);
    return ( (std::uint64_t)a << 32 ) | (std::uint64_t)b;
}


// euclidean distance between two cells using their positions
static inline double dist_cells(Cell* A, Cell* B)
{
    double dx = A->position[0] - B->position[0];
    double dy = A->position[1] - B->position[1];

    double dz = A->position[2] - B->position[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// attach two cells and record their current separation as their rest length 
static inline void attach_cells_with_rest(Cell* A, Cell* B)
{
    attach_cells(A,B);

    std::uint32_t ida = (std::uint32_t) A->ID;
    std::uint32_t idb = (std::uint32_t) B->ID;

    g_rest_length[pair_key(ida,idb)] = dist_cells(A,B);
}


// look up rest length for two cells and return fall back value if it doesn't exist
static inline double get_rest_length(Cell* A, Cell* B, double fallback)
{
    std::uint32_t ida = (std::uint32_t) A->ID;

    std::uint32_t idb = (std::uint32_t) B->ID;

    auto it = g_rest_length.find(pair_key(ida,idb));
    if(it != g_rest_length.end())
        return it->second;

    return fallback;
}



// Helper functions to ensure initial rod placement doesn't intersect or overlap for 2d

struct RodSegment { double x1,y1,x2,y2; };

static double orient2d(double ax,double ay,double bx,double by,double cx,double cy)
{
    // cross product between AB and AC tells us which side of line AB
    // AC lies on 

    // if result > 0, point C lies to the left 
    // if result < 0, point C lies to the right 
    // if result = 0, A, B and C lie on the same line 
    return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
}

static bool on_segment2d(double ax,double ay,double bx,double by,double px,double py)
{
    // Assuming P lies on the same line segment as AB, 
    // this checks whether P lies in between A and B 
    // or beyond it
    return (std::min(ax,bx) <= px && px <= std::max(ax,bx) &&
            std::min(ay,by) <= py && py <= std::max(ay,by));
}

// True intersection test (catches X crossings)
static bool segments_intersect2d(double ax,double ay,double bx,double by,
                                 double cx,double cy,double dx,double dy)
{
    // which side of AB is point C
    double o1 = orient2d(ax,ay,bx,by,cx,cy);
    // which side of AB is point B

    // Note: two segements intersect if C is on one side of ABD and B is 
    // on the other side of AB
    // and the same goes for the other line segment 

    // as such o1 and o2 should be different
    // and o3 and o4 should also be different 
    double o2 = orient2d(ax,ay,bx,by,dx,dy);
    double o3 = orient2d(cx,cy,dx,dy,ax,ay);
    double o4 = orient2d(cx,cy,dx,dy,bx,by);

    // if they are different, it shows, there is no intersection
    if( (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0) ) return true;


    // a point can be colinear but if it outside the segment
    // then there is no intersection 
    if( o1 == 0 && on_segment2d(ax,ay,bx,by,cx,cy) ) return true;
    if( o2 == 0 && on_segment2d(ax,ay,bx,by,dx,dy) ) return true;
    if( o3 == 0 && on_segment2d(cx,cy,dx,dy,ax,ay) ) return true;
    if( o4 == 0 && on_segment2d(cx,cy,dx,dy,bx,by) ) return true;

    return false;
}





static double dot2(double ax,double ay,double bx,double by)
{
    // dot product between two vectors 
    // helps us figure out how far along a point the nearest point is
    return ax*bx + ay*by;
}

static double dist2_point_to_segment(double px,double py,
                                     double ax,double ay,double bx,double by)
{

    // given point P (px, py) and segment end points 
    // A(ax, ay), B(bx, by), find the closest point on the segment to P
    // return the squared distance 
    double abx = bx-ax, aby = by-ay; // AB vector
    double apx = px-ax, apy = py-ay; // AP vector
    double ab2 = abx*abx + aby*aby; // |AB}^2 to normalise projection 

    // t tells is how long along AB the closest point is
    // t = 0, closest point is A 
    // t = 1, closest point is B 
    double t = 0.0;

    // project AP into AB to find closest point 
    if(ab2 > 0.0) t = dot2(apx,apy,abx,aby)/ab2;

    // ensure t lies within 0 and 1
    t = std::max(0.0, std::min(1.0, t));

    // compute actual closest point
    double cx = ax + t*abx;
    double cy = ay + t*aby;

    // find vector from closest point C to P 
    double dx = px - cx;
    double dy = py - cy;

    // return the squared distance
    return dx*dx + dy*dy;
}

// Safe distance check once intersection is excluded
static double dist2_segment_to_segment(double ax,double ay,double bx,double by,
                                       double cx,double cy,double dx,double dy)
{

    // if two lines segments don't cross, then the closest points must involve alteast on end point 

    // distance from A and B to segment CD
    double d1 = dist2_point_to_segment(ax,ay,cx,cy,dx,dy);
    double d2 = dist2_point_to_segment(bx,by,cx,cy,dx,dy);

    // distance from C and D to segment AB
    double d3 = dist2_point_to_segment(cx,cy,ax,ay,bx,by);
    double d4 = dist2_point_to_segment(dx,dy,ax,ay,bx,by);

    // return the smallest of the four end point to segement distances 
    return std::min(std::min(d1,d2), std::min(d3,d4));
}



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
    cell_defaults.functions.custom_cell_rule      = NULL;

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
        pT->functions.custom_cell_rule = avoid_boundaries;
    }

    // rod cell
    Cell_Definition* pRod = find_cell_definition( "rod" );
    if( pRod == nullptr )
    {
        std::cout << "no rod definition\n";
    }
    else
    {
        pRod->phenotype.motility.is_motile = false;

        // raise attachment limit
        int idx = pRod->custom_data.find_variable_index( "max_attachments" );
        if( idx < 0 )
            pRod->custom_data.add_variable( "max_attachments", "dimensionless", 300);
        else
            pRod->custom_data[ "max_attachments" ] = 300;

        // contact spring for rod rod attachments
        pRod->functions.contact_function = contact_function;

        // IL-2 secretion
        pRod->functions.update_phenotype = secretion_rate_rod_cells;
    }

    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    setup_cell_rules();
    display_cell_definitions( std::cout );
}


void setup_microenvironment( void )
{
    initialize_microenvironment();
}


void setup_tissue( void )
{
    // Domain
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];

    if( default_microenvironment_options.simulate_2D )
    {
        Zmin = 0.0;
        Zmax = 0.0;
    }

    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;

    load_cells_from_pugixml();
    g_rest_length.clear();

    // T cells 
    Cell_Definition* pT = find_cell_definition("T cell");
    if(pT)
    {
        // int N = compute_number_of_tcells_from_density();
        int N = parameters.ints("number_of_T_cells");

        for(int i=0; i<N; i++)
        {
            Cell* c = create_cell(*pT);
            double x = Xmin + UniformRandom()*Xrange;
            double y = Ymin + UniformRandom()*Yrange;
            double z = default_microenvironment_options.simulate_2D ? 0.0
                      : Zmin + UniformRandom()*Zrange;
            c->assign_position(x,y,z);
        }
    }

    // rod parameters 
    // Rods are arranged according to thin hexagonal close packed model
    // Each rod is a stack of "layers" along the rod axis.
    // Each layer contains exactly 3 cells forming a triangle
    // Then we ABAB-stagger layers where:
    // layer A: triangle at offsets t0,t1,t2
    // layer B: same triangle shifted by shiftB so cells sit in the gaps
    // the two layers are repeated along the rod axis
 

    Cell_Definition* pRodDef = find_cell_definition("rod");
    if(!pRodDef) return;

    int number_of_rods = compute_number_of_rods_from_concentration();

    if(number_of_rods <= 0)
        number_of_rods = parameters.ints("number_of_rods");

    int cells_per_rod = parameters.ints("rod_cells_per_rod");

    // Spacing between rod cells 
    double rad = pRodDef->phenotype.geometry.radius;
    // d is the in layer nearest neighbour spacing 
    // it is consistent with the built in PhysiCell confluent length scaling
    double d   = 2.0 * rad * 0.95238095238;  
    // ds is the vertical spacing between consecutive layers along the rod axis
    // use pythag
    double ds  = std::sqrt(2.0/3.0) * d;          

    int layers_per_rod = std::max(2, cells_per_rod / 3);
    double L = (layers_per_rod - 1) * ds; // rod length 
    double halfL = 0.5 * L;

    // rod placement
    double minDist  = 3.0 * d;
    double minDist2 = minDist * minDist;
    double margin   = minDist;
    int max_tries   = 2000;


    // this contains all the rods
    // std::vector<RodSegment> placed_rods;
    // placed_rods.reserve(number_of_rods);

    // Determine the z values for the rod layers based on the number of layer 

    const int NUM_LEVELS = 3; //TODO: make parameter in xml

    std::vector<RodSegment> placed_rods_layer[NUM_LEVELS];
    for(int L = 0; L < NUM_LEVELS; L++)
    {
        placed_rods_layer[L].reserve(number_of_rods);
    }

    double z_layer_sep = (std::sqrt(3.0) / 2.0) * d + 7.0 * rad; // distance between rod layers 
    double z_levels[NUM_LEVELS];

    // if 2d, z = 0
    if(default_microenvironment_options.simulate_2D)
    {
        for(int L = 0; L < NUM_LEVELS; L++)
            z_levels[L] = 0.0;
    }
    else
    {
        // compute z values for layers 
        for(int L = 0; L < NUM_LEVELS; L++)
        {
            z_levels[L] = (L - 0.5 * (NUM_LEVELS - 1)) * z_layer_sep;
        }


        for(int L = 0; L < NUM_LEVELS; L++)
        {
            if(z_levels[L] < Zmin + 0.5 * z_layer_sep)
                z_levels[L] = Zmin + 0.5 * z_layer_sep;

            if(z_levels[L] > Zmax - 0.5 * z_layer_sep)
                z_levels[L] = Zmax - 0.5 * z_layer_sep;
        }
    }

    std::cout << "z_levels: ";
    for(int L = 0; L < NUM_LEVELS; L++)
        std::cout << z_levels[L] << " ";
    std::cout << std::endl;


    int rods_created = 0;

    for(int rod_idx=0; rod_idx<number_of_rods; rod_idx++)
    {
        bool placed = false;

        double x0,y0,z0,theta,ux,uy,vx,vy;
        double ax,ay,bx,by;

        int chosen_level = -1;

        for(int attempt=0; attempt<max_tries; attempt++)
        {
            // random centre (x,y only)
            x0 = Xmin + UniformRandom()*Xrange;
            y0 = Ymin + UniformRandom()*Yrange;

            // random orientation in the x-y plane
            theta = 2.0*M_PI*UniformRandom();
            ux = cos(theta); uy = sin(theta);
            vx = -uy;        vy = ux;

            // rod end points in XY
            ax = x0 - halfL*ux; ay = y0 - halfL*uy;
            bx = x0 + halfL*ux; by = y0 + halfL*uy;

            // keep within XY boundaries
            if(ax < Xmin+margin || ax > Xmax-margin ||
            bx < Xmin+margin || bx > Xmax-margin ||
            ay < Ymin+margin || ay > Ymax-margin ||
            by < Ymin+margin || by > Ymax-margin)
                continue;

            // try place in z level 1, if it intersecs then level 2 etc...
            for(int L=0; L<NUM_LEVELS; L++)
            {
                bool ok = true;

                for(const auto& r : placed_rods_layer[L])
                {
                    if( segments_intersect2d(ax,ay,bx,by, r.x1,r.y1,r.x2,r.y2) )
                    { ok = false; break; }

                    double d2 = dist2_segment_to_segment(ax,ay,bx,by, r.x1,r.y1,r.x2,r.y2);
                    if(d2 < minDist2)
                    { ok = false; break; }
                }

                if(ok)
                {
                    placed = true;
                    chosen_level = L;

                    // FIXED z for that layer
                    z0 = z_levels[L];

                    // record this rod only in that layer
                    placed_rods_layer[L].push_back({ax,ay,bx,by});
                    break;
                }
            }

            if(placed) break;
        }




        if(!placed)
        {
            std::cout << "could not place rod " << rod_idx << "\n";
            continue;
        }

        // u points towards the rod 
        Vec3 u = normalize3({ux,uy,0.0});
        // v lies orthogonal to the rod
        Vec3 v = normalize3({vx,vy,0.0});
        // w lies straights up
        Vec3 w = {0.0,0.0,1.0};

        // centre 
        Vec3 center = {x0,y0,z0};


        // t0, t1 and t2 represent an equilaterial triangle offsets of length d in the v-w plane from the centre 
        Vec3 t0 = {0,0,0};
        Vec3 t1 = d * v;
        Vec3 t2 = (0.5*d) * v + (std::sqrt(3.0)/2.0 * d) * w;

        // (A + B + C ) / 3 is the centroid of a triangle 
        Vec3 shiftB = (1.0/3.0) * (t1 + t2); 

        std::vector<Cell*> rod_cells;
        rod_cells.reserve(3 * layers_per_rod);


        for(int k=0; k<layers_per_rod; k++)
        {
            double s = -halfL + k*ds; // ds is the vertical spacing and -halfL half the length of the rod

            // if k is off apply a shift, otherwise don't 
            Vec3 shift = (k % 2 == 0) ? Vec3{0,0,0} : shiftB;
            Vec3 offs[3] = { t0 + shift, t1 + shift, t2 + shift };

            for(int m=0; m<3; m++)
            {
                // from the centre of the rod, move s away in the direction of u 
                // offset is applied every second layer since the cells fill the gaps of the previous layer 
                Vec3 P = center + s*u + offs[m]; 

                // assign the x, y and z  position of the rod cells 
                Cell* c = create_cell(*pRodDef);
                c->assign_position(P.x, P.y, P.z);

                c->custom_data["rod_id"] = rod_idx;
                c->custom_data["row"]    = k;
                c->custom_data["col"]    = m;

                rod_cells.push_back(c);
            }
        }

        // Attachments:
        // triangle edges in each layer 
        for(int k=0; k<layers_per_rod; k++)
        {
            Cell* a = rod_cells[3*k + 0];
            Cell* b = rod_cells[3*k + 1];
            Cell* c = rod_cells[3*k + 2];

            attach_cells_with_rest(a,b);
            attach_cells_with_rest(b,c);
            attach_cells_with_rest(c,a);
        }

        // nearest neighbours between adjacent layers (distance close to d)
        const double tol = 0.08 * d;
        for(int k=0; k<layers_per_rod-1; k++)
        {
            for(int m=0; m<3; m++)
            {
                Cell* c0 = rod_cells[3*k + m];
                for(int n=0; n<3; n++)
                {
                    Cell* c1 = rod_cells[3*(k+1) + n];

                    double dxp = c0->position[0] - c1->position[0];
                    double dyp = c0->position[1] - c1->position[1];
                    double dzp = c0->position[2] - c1->position[2];
                    double dist = std::sqrt(dxp*dxp + dyp*dyp + dzp*dzp);

                    if(std::fabs(dist - d) <= tol)
                        attach_cells_with_rest(c0, c1);
                }
            }
        }




        // longitudinal braces 
        for(int k=0; k<layers_per_rod-2; k++)
        {
            for(int m=0; m<3; m++)
            {
                Cell* c0 = rod_cells[3*k + m];
                Cell* c2 = rod_cells[3*(k+2) + m];
                attach_cells_with_rest(c0, c2);
            }
        }

        // // diagonal braces  (keep it off for now)
        // // so for each cell in k, attach to the cell in k+1 that isn't of the same index
        // for(int k=0; k<layers_per_rod-1; k++)
        // {
        //     for(int m=0; m<3; m++)
        //     {
        //         Cell* c0 = rod_cells[3*k + m];

        //         Cell* cA = rod_cells[3*(k+1) + ((m+1)%3)];
        //         Cell* cB = rod_cells[3*(k+1) + ((m+2)%3)];

        //         attach_cells_with_rest(c0, cA);
        //         attach_cells_with_rest(c0, cB);
        //     }
        // }

        rods_created++;
    }

    std::cout << "Rods requested=" << number_of_rods
              << " rods created=" << rods_created << std::endl;

    for(int L=0; L<NUM_LEVELS; L++)
    std::cout << "Layer " << L << " rods: " << placed_rods_layer[L].size() << "\n";

}


static void rod_elastic_contact_thin_hcp(
    Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
    if( pC1->position.size() != 3 || pC2->position.size() != 3 ) return;

    static int rod_type = get_cell_definition("rod").type;
    if(pC1->type != rod_type || pC2->type != rod_type) return;

    if(!are_attached(pC1, pC2)) return;

    std::vector<double> displacement = pC2->position;
    displacement -= pC1->position;

    int ii = find_cell_definition_index( pC1->type );
    int jj = find_cell_definition_index( pC2->type );

    double adhesion_ii = pC1->phenotype.mechanics.attachment_elastic_constant *
                         pC1->phenotype.mechanics.cell_adhesion_affinities[jj];
    double adhesion_jj = pC2->phenotype.mechanics.attachment_elastic_constant *
                         pC2->phenotype.mechanics.cell_adhesion_affinities[ii];

    double effective_k = std::sqrt( adhesion_ii * adhesion_jj );

    // fallback rest length
    double d0 = 2.0 * p1.geometry.radius * 0.95238095238;

    // per edge rest length 
    double rest = get_rest_length(pC1, pC2, d0);

    double strength = ( norm(displacement) - rest ) * effective_k;
    normalize( &displacement );
    axpy( &(pC1->velocity), strength, displacement );
}

void contact_function( Cell* pMe, Phenotype& phenoMe,
                       Cell* pOther, Phenotype& phenoOther, double dt )
{
    rod_elastic_contact_thin_hcp(pMe, phenoMe, pOther, phenoOther, dt);
}


void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
    return;
}

void tcell_update_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    standard_update_cell_velocity( pCell, phenotype, dt );

    static int rod_type = get_cell_definition("rod").type;

    int k_touch = pCell->custom_data.find_variable_index("touching_rod");
    int k_time  = pCell->custom_data.find_variable_index("attached_time");
    int k_state = pCell->custom_data.find_variable_index("state");

    if( k_touch < 0 || k_time < 0 || k_state < 0 ) return;

    pCell->custom_data[k_touch] = 0.0;

    auto nearby = pCell->nearby_interacting_cells();

    for( Cell* other : nearby )
    {
        if( other->type != rod_type ) continue;

        double dx = pCell->position[0] - other->position[0];
        double dy = pCell->position[1] - other->position[1];
        double dz = pCell->position[2] - other->position[2];
        double d  = std::sqrt(dx*dx + dy*dy + dz*dz);

        double contact_dist = pCell->phenotype.geometry.radius + other->phenotype.geometry.radius;

        if( d <= contact_dist )
        {
            pCell->custom_data[k_touch] = 1.0;
            break;
        }
    }

    if( pCell->custom_data[k_state] < 0.5 )
    {
        if( pCell->custom_data[k_touch] > 0.5 )
        {
            pCell->custom_data[k_time] += dt;

            const double THRESH_MIN = 30;

            if( pCell->custom_data[k_time] >= THRESH_MIN )
            {
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
    if( pParent == nullptr || pDaughter == nullptr ) return;

    static int t_type = get_cell_definition("T cell").type;
    if( pParent->type != t_type ) return;

    int k_state = pParent->custom_data.find_variable_index("state");
    int k_time  = pParent->custom_data.find_variable_index("attached_time");
    int k_touch = pParent->custom_data.find_variable_index("touching_rod");

    if( k_state < 0 ) return;

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

    clamp_cell_to_domain(pParent);
    clamp_cell_to_domain(pDaughter);

}


std::vector<std::string> coloring_function( Cell* pCell )
{
    std::vector<std::string> output(4, "grey");

    if( pCell->phenotype.death.dead )
    {
        output[0] = "red";
        output[2] = "darkred";
        return output;
    }

    static int t_type   = get_cell_definition("T cell").type;
    static int rod_type = get_cell_definition("rod").type;

    if( pCell->type == t_type )
    {
        int k_state = pCell->custom_data.find_variable_index("state");
        bool active = (k_state >= 0 && pCell->custom_data[k_state] > 0.5);

        if( active )
        {
            output[0] = "#FF8C00";
            output[2] = "#FFFFFF";
        }
        else
        {
            output[0] = "blue";
            output[2] = "darkblue";
        }
    }
    else if( pCell->type == rod_type )
    {
        output[0] = "green";
        output[2] = "darkgreen";
    }
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
            static int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
            static int cycle_end_index   = live.find_phase_index(PhysiCell_constants::live);

            static int IL2_index = microenvironment.find_density_index("IL-2");
            double IL2 = pCell->nearest_density_vector()[IL2_index];

            double rPmax = parameters.doubles("rPmax");
            double IP    = parameters.doubles("IP");

            phenotype.cycle.data.transition_rate( cycle_start_index, cycle_end_index ) = rPmax*IL2/(IP+IL2);
        }
    }
}

void secretion_rate_rod_cells( Cell* pCell , Phenotype& phenotype, double dt )
{
    static int IL2_index = microenvironment.find_density_index("IL-2");

    double v   = parameters.doubles("v");
    double q   = parameters.doubles("q");
    double rho = parameters.doubles("rho");

    double time = PhysiCell::PhysiCell_globals.current_time;

    phenotype.secretion.secretion_rates[IL2_index] = v*q*rho*exp(-v*time);
}


int compute_number_of_rods_from_concentration()
{
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
    double domain_volume_mL  = domain_volume_um3 * 1e-12;

    double total_mass_ug = conc_ug_per_mL * domain_volume_mL;

    Cell_Definition* pRodDef = find_cell_definition("rod");
    double rod_cell_vol_um3  = pRodDef->phenotype.volume.total;

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



void avoid_boundaries( Cell* pCell )
{
	// add velocity to steer clear of the boundaries 
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static double avoid_zone = 1; 
	static double avoid_speed = -0.05; 
	
	// near edge: 
	bool near_edge = false; 
	if( pCell->position[0] < Xmin + avoid_zone || pCell->position[0] > Xmax - avoid_zone )
	{ near_edge = true; } 
	
	if( pCell->position[1] < Ymin + avoid_zone || pCell->position[1] > Ymax - avoid_zone )
	{ near_edge = true; } 
	
	if( default_microenvironment_options.simulate_2D == false )
	{
		if( pCell->position[2] < Zmin + avoid_zone || pCell->position[2] > Zmax - avoid_zone )
		{ near_edge = true; } 
	}
	
	if( near_edge )
	{
		pCell->velocity = pCell->position; // move towards origin 
		pCell->velocity *= avoid_speed; // move towards origin 
	}
	
	return; 
}

void avoid_boundaries( Cell* pCell , Phenotype& phenotype, double dt )
{ return avoid_boundaries( pCell ); } 




void clamp_cell_to_domain(Cell* c)
{
    double domain_x_min = microenvironment.mesh.bounding_box[0];
    double domain_y_min = microenvironment.mesh.bounding_box[1];
    double domain_z_min = microenvironment.mesh.bounding_box[2];

    double domain_x_max = microenvironment.mesh.bounding_box[3];
    double domain_y_max = microenvironment.mesh.bounding_box[4];
    double domain_z_max = microenvironment.mesh.bounding_box[5];

    if (c->position[0] < domain_x_min) c->position[0] = domain_x_min;
    if (c->position[0] > domain_x_max) c->position[0] = domain_x_max;

    if (c->position[1] < domain_y_min) c->position[1] = domain_y_min;
    if (c->position[1] > domain_y_max) c->position[1] = domain_y_max;

    if (c->position[2] < domain_z_min) c->position[2] = domain_z_min;
    if (c->position[2] > domain_z_max) c->position[2] = domain_z_max;
}




















// void clamp_cell_to_domain(Cell* c)
// {

//     auto& bb = microenvironment.mesh.bounding_box;

//     c->position[0] = std::max(bb[0], std::min(bb[3],  c->position[0]));
//     c->position[1] = std::max(bb[1], std::min(bb[4], c->position[1]));
//     c->position[2] = std::max(bb[2], std::min(bb[5], c->position[2]));

// }

