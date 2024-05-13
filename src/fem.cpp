#include "fem.h"
#include "mesh.h"
#include "simu.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
    	//std::cout << "Coordonnee des sommets : \n";
        if ( border==true ) {
        	for (int v = 0; v < 2; ++v ){
        		vertices_.push_back(M.get_edge_vertex(i,v));
        		//std::cout << "v" << v <<" : " << vertices_[v].x << "	" << vertices_[v].y << std::endl;
        	} 
        }
        else {
        	for ( int v = 0; v < 3; ++v ){
        		vertices_.push_back(M.get_triangle_vertex(i,v)) ;
        		//std::cout << "v" << v <<" : " << vertices_[v].x << "	" << vertices_[v].y << std::endl;
        	} 
        }
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        vertex r ;
        if ( border_ ) {
        	double xi = x_r.x ;
        	r.x = (1-xi)*vertices_[0].x + xi*vertices_[1].x ;
        	r.y = (1-xi)*vertices_[0].y + xi*vertices_[1].y ;
        }
        else {
        	double xi = x_r.x ;
        	double eta = x_r.y ;
        	r.x = (1-xi-eta)*vertices_[0].x + xi*vertices_[1].x + eta*vertices_[2].x ;
        	r.y = (1-xi-eta)*vertices_[0].y + xi*vertices_[1].y + eta*vertices_[2].y ;
		}
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        DenseMatrix J ;
        if ( border_ ) {
        	J.set_size(2,1);
        	J.set(0,0, -vertices_[0].x + vertices_[1].x) ;
        	J.set(1,0, -vertices_[0].y + vertices_[1].y) ;
        }
        else {
        	J.set_size(2,2);
        	J.set(0,0, -vertices_[0].x + vertices_[1].x) ;
        	J.set(0,1, -vertices_[0].x + vertices_[2].x) ;
        	J.set(1,0, -vertices_[0].y + vertices_[1].y) ;
        	J.set(1,1, -vertices_[0].y + vertices_[2].y) ;
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {

	DenseMatrix J = jacobian_matrix(x_r);
	if (border_) {
		double JTJ = J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0);
		return std::sqrt(JTJ);
	}
        else {
        	double JTJ = J.det_2x2();
        	return JTJ;
        }
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        bool SF_constru = true;
        if (dim != 1 && dim != 2) {
        	std::cout << "Implementation uniquement en 1D ou 2D \n";
        	SF_constru = false;
        } 
        if (order != 1) {
        	std::cout << "Uniquement fonctions d'ordre 1 \n";
        	SF_constru = false;
        }
        assert(SF_constru);
    }

    int ShapeFunctions::nb_functions() const
    {
        if (dim_ == 1) {
        	return 2 ;
        }
        if (dim_ == 2) {
        	return 3;
        }
        return 0;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {	
    	double g;
    	if (dim_ == 1 ) {
		double xi = x_r.x;
		switch(i) {
			case (0):
				g = 1-xi; break;
			case (1):
				g = xi; break;
		}
    	} 
    	else {
		double xi = x_r.x;
		double eta = x_r.y;
		switch(i) {
			case (0):
				g = 1-xi-eta; break;
			case (1):
				g = xi; break;
			case (2):
				g = eta; break;
		}
    	}
    	return g;
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        vec2 g ;
        if (dim_ == 1 ) {
        	switch(i) {
        		case (0):
        			g.x = -1.; break;
        		case (1):
        			g.x = 1.; break;
        	}
        	g.y = 0. ;
        } 
        else {
        	switch(i) {
        		case (0):
        			g.x = -1.; g.y = -1.; break;
        		case (1):
        			g.x = 1.; g.y = 0.; break;
        		case (2):
        			g.x = 0.; g.y = 1.; break;
        	}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex), 
        DenseMatrix& Ke )
    {
    	Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());
    	for(int i = 0; i<reference_functions.nb_functions();++i){
    		for (int j = 0; j<reference_functions.nb_functions();++j){
    			double somme = 0;
    			for (int q = 0; q < quadrature.nb_points(); ++q){
    				vertex p_q = quadrature.point(q);
    				double w_q = quadrature.weight(q);
    				DenseMatrix J_inv = elt_mapping.jacobian_matrix(p_q).invert_2x2();
    				vec2 grad_i = J_inv.transpose().mult_2x2_2(reference_functions.evaluate_grad(i, p_q));
    				vec2 grad_j = J_inv.transpose().mult_2x2_2(reference_functions.evaluate_grad(j, p_q));
    				somme += w_q * coefficient(elt_mapping.transform(p_q)) * dot(grad_i, grad_j) * elt_mapping.jacobian(p_q);
    			}
    			Ke.set(i,j, somme);
    		}
    	}
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
    	for (int i = 0; i<Ke.height(); ++i) {
    		for (int j = 0; j<Ke.width(); ++j){
    			int I = M.get_triangle_vertex_index(t, i);
    			int J = M.get_triangle_vertex_index(t, j);
    			K.add(I,J, Ke.get(i,j));
    		}
    	}
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
    	for(int i = 0; i<reference_functions.nb_functions();++i){
    		double somme = 0;
    		for (int q = 0; q < quadrature.nb_points(); ++q){
    			vertex p_q = quadrature.point(q);
    			double w_q = quadrature.weight(q);
    			somme += w_q * source(elt_mapping.transform(p_q)) * reference_functions.evaluate(i, p_q) * elt_mapping.jacobian(p_q);
    		} 
    		Fe[i] = somme;
    	}
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
    	for(int i = 0; i<reference_functions_1D.nb_functions();++i){
    		double somme = 0;
    		for (int q = 0; q < quadrature_1D.nb_points(); ++q){
    			vertex p_q = quadrature_1D.point(q);
    			double w_q = quadrature_1D.weight(q);
    			somme += w_q * neumann(elt_mapping_1D.transform(p_q)) * reference_functions_1D.evaluate(i, p_q) * elt_mapping_1D.jacobian(p_q);
    		} 
    		Fe[i] = somme;
    	}
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
    	for (int j = 0; j<Fe.size(); ++j){
    		int J = M.get_triangle_vertex_index(i, j);
    		F[J] += Fe[j];
    		}
    }
    
    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::vector< bool > processed_vertices(values.size(), false);
        double penalty = 10000.;
        for( int edge = 0; edge < M.nb_edges(); edge++ ) {
        	int edge_attribute = M.get_edge_attribute(edge);
        	if( attribute_is_dirichlet[edge_attribute] ) {
        	        for( int v = 0; v < 2; v++ ) {
        	        	int vertex_index = M.get_edge_vertex_index(edge, v);
        	        	if( !processed_vertices[vertex_index] ) {
        	                	processed_vertices[vertex_index] = true;
                        		K.add(vertex_index, vertex_index, penalty);
                        		F[vertex_index] += penalty*values[vertex_index];
                    }
                }
            }
        }
    }

    void solve_poisson_problem(
            const Mesh* M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {        
    std::cout << "solve poisson problem" << '\n';
    }

}
