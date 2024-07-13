#include "linesearch_rcg.h"

namespace ModuleDirectMin
{

    void LineSearchRCG::set_default_parameters()
    {
        LineSearchBase::set_default_parameters();

        sigma = 0.0;

        cg_algo = DAI_YUAN;

        cg_algo_names = new std::string[RCG_ALGO_LENGTH];

        cg_algo_names[FLETCHER_REEVES].assign("FLETCHER_REEVES");
        cg_algo_names[DAI_YUAN].assign("DAI_YUAN");
        cg_algo_names[POLAK_RIBIERES].assign("POLAK_RIBIERES");
        cg_algo_names[POLAK_RIBIERES_MOD].assign("POLAK_RIBIERES_MOD");
        cg_algo_names[FR_PR].assign("FR_PR");
        cg_algo_names[HAGER_ZHANG].assign("HAGER_ZHANG");
        cg_algo_names[HESTENES_STIEFEL].assign("HESTENES_STIEFEL");

        method_name.assign("RCG"); 
    }

    void LineSearchRCG::update_parameters(LineSearchOptions* ls_opt_in)
    {
        LineSearchBase::update_parameters(ls_opt_in);
        set_cg_algo(ls_opt_in->ls_cg_algo); // set the cg algorithm using input option

    }

    void LineSearchRCG::get_search_direction()
    {

        // FLT_EPSILON      = 1.19209e-07
        // DBL_EPSILON      = 2.22045e-16
        // LDBL_EPSILON     = 1.0842e-19
        if( iter == 0 || 
            x1.inner_product(eta1, gf1) / ngf1 /ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon() )
          )
            // metric(x1, eta1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon() )

        {
            Pgf1 = prob->preconditioner(gf1); //Use preconditioner
            eta1 = -Pgf1;

            // if(iter > 0 && verbose)
            // {
            //     // eta1.brief_print("inside get_search_direction, eta1:");
            //     std::cout << "metric(x1, eta1, gf1) / ngf1 / ngf1 = " 
            //              << metric(x1, eta1, gf1) / ngf1 / ngf1 
            //              << std::endl;
            // }

            if( x1.inner_product(eta1, gf1) / ngf1 /ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon() )
            )
// metric(x1, eta1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon() ))
            {   
                eta1 = -1.0 * gf1; // Without preconditioner
            }
        }
    }


    void LineSearchRCG::update_data()
    {
        // compute_sigma();
        // Stiefel Td1, Tgf1;

        // Td1 = vectran(x2, eta1);
        double sigmaFR, sigmaPR;

        Stiefel Td1, Tgf1;
        Stiefel y1;
        Td1 = vectran(x2, eta1);  nV ++;
        Tgf1 = vectran(x2, gf1); nV ++;

        Pgf2 = prob->preconditioner(gf2);
        // Stiefel Td1, Tgf1;
        // Stiefel y1;

        if(fabs(metric(x2, gf2, Tgf1)) / metric( x1, gf1, gf1) > 0.1)
        { 
            /*Restart using the safeguard (5.52) in [NW06] */
            sigma = 0;
        }
        else
        {
            switch(cg_algo)
            {
                case FLETCHER_REEVES:
                    sigma = metric(x2, gf2, Pgf2) 
                                /metric(x1, gf1, Pgf1);
                    break;
                case DAI_YUAN:
                    y1 = gf2 - Tgf1;
                    sigma = metric(x1, gf2, Pgf2) 
                                /metric(x2, Td1, y1);
                    break;
                case POLAK_RIBIERES:
                    y1 = gf2 - Tgf1;
                    sigma = metric(x2, Pgf2, y1) 
                                /metric(x1, gf1, Pgf1);
                    break;

                // case POLAK_RIBIERES_MOD:
                //     y1 = gf2 - Tgf1;
                //     sigma = metric(x2, y1, gf2) 
                //             /metric(x1, gf1, gf1);
                //     sigma = (sigma < 0) ? 0 : sigma;
                //     break;

                case FR_PR:
                    sigmaFR = metric(x2, gf2, gf2) 
                                /metric(x1, gf1, gf1);
                    sigmaPR = metric(x2, gf2, Tgf1) 
                                /metric(x1, gf1, gf1);
                    sigma = (sigmaPR < -sigmaFR) ? - sigmaFR : ((sigmaPR > sigmaFR) ? sigmaFR : sigmaPR);
                    break;

                // case HAGER_ZHANG:
                //     y1 = gf2 - Tgf1;
                //     sigma = metric(x2, Tgf1, Pgf2) / metric 
                //     break;

                case HESTENES_STIEFEL:
                    y1 = gf2 - Tgf1;
                    sigma = metric(x2, Pgf2, Tgf1) 
                                /metric(x2, Tgf1, Td1);
                    break;
                default:
                    std::cout << "Unknown RCG algorithm" << std::endl;
                    sigma = 0;
            }
                    
        }
        // if (cg_algo == FLETCHER_REEVES)
        // {
        //     Td1 = vectran(x2, eta1);
        //     sigma = metric(x2, gf2, gf2) 
        //                 /metric(x1, gf1, gf1);
        // }
        // else if (cg_algo == POLAK_RIBIERES ) 
        // {
        //     Td1 = vectran(x2, eta1);
        //     sigma = metric(x2, Td1, gf2) 
        //                 /metric(x1, gf1, gf1);
        // }
        // else if (cg_algo == HESTENES_STIEFEL)
        // {
        //     Td1 = vectran(x2, eta1);
        //     Tgf1 = vectran(x2, gf1);

        //     y1 = gf2 - Tgf1;
        //     sigma = metric(x1, gf2, y1) 
        //                 /metric(x2, Td1, y1);
        // }

        // else if (cg_algo == DAI_YUAN)
        // {

        //     Td1 = vectran(x2, eta1);
        //     Tgf1 = vectran(x2, gf1);
        //     y1 = gf2 - Tgf1;
        //     sigma = metric(x1, gf2, gf2) 
        //                 /metric(x2, Td1, y1);
        // }

        Pgf1 = Pgf2;
        if (sigma == 0)
        {
            eta1 = -gf2;
        }
        else 
        {
            eta1 = - gf2 + sigma * Td1; // update new direction
        }
        // std::cout << "sigma=" << sigma << std::endl;
    }


    void LineSearchRCG::set_cg_algo( std::string & name)
    {
        if (name == "fr")
        {
            cg_algo = FLETCHER_REEVES;
        }
        else if (name == "dy")
        {
            cg_algo = DAI_YUAN;
        }
        else if (name == "pr")
        {
            cg_algo = POLAK_RIBIERES;
        }  
        else if (name == "pr_mod")
        {
            cg_algo = POLAK_RIBIERES_MOD;
        }    
        else if (name == "hz")
        {
            cg_algo = HAGER_ZHANG;
        }
        else if (name == "hs")
        {
            cg_algo = HESTENES_STIEFEL;
        }
        else if (name == "fr_pr")
        {
            cg_algo = FR_PR;
        }   
        else
        {
            std::cerr << "Unknown CG algorithm: " << name << std::endl;
            exit(1);
        }
    }
    
    void LineSearchRCG::compute_sigma()
    {
        Stiefel Td1, Tgf1;
        Stiefel y1;
        Td1 = vectran(x2, eta1);  nV ++;
        Tgf1 = vectran(x2, gf1); nV ++;

        if(fabs(metric(x2, gf2, Tgf1)) / metric( x1, gf1, gf1) > 0.1)
        {
            sigma = 0;
        }
        switch(cg_algo)
        {
            case FLETCHER_REEVES:
                sigma = metric(x2, gf2, gf2) 
                            /metric(x1, gf1, gf1);
                break;
            case DAI_YUAN:
                y1 = gf2 - Tgf1;
                sigma = metric(x1, gf2, gf2) 
                            /metric(x2, Td1, y1);
                break;
            case POLAK_RIBIERES:
                sigma = metric(x2, Td1, gf2) 
                        /metric(x1, gf1, gf1);
                break;
            case HAGER_ZHANG:
                break;
            case HESTENES_STIEFEL:
                y1 = gf2 - Tgf1;
                sigma = metric(x1, Pgf2, y1) 
                            /metric(x2, Td1, y1);
                break;
            default:
                std::cout << "Unknown RCG algorithm" << std::endl;
                sigma = 0;
        }
    }


} // namespace ModuleDirectMin