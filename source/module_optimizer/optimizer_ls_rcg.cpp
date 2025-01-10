#include "optimizer_ls_rcg.h"

namespace Module_Optimizer
{
    // RCG::RCG()
    // {

    //     // OptimizerLSBase();
    //     SetDefaultParams();

    // }
    // RCG::RCG(Problem * prob)
    // {
    //     std::cout << "Inside RCG::RCG" << std::endl;
    //     SetDefaultParams();
    //     initialize(prob);
    //     // set_cg_algo(algoName);
    // }

    RCG::RCG(Problem *prob_in, LSOptions * opt_in)
    {
        set_default_params();
        // OptimizerLSBase::initialize(prob_in);
        OptimizerLSBase::initialize();

        update_params(opt_in);
    }

    void RCG::set_default_params()
    {
        OptimizerLSBase::set_default_params();
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

        name.assign("RCG"); 

    }
    void RCG::update_params(LSOptions* opt_in)
    {
        OptimizerLSBase::update_params(opt_in);
        set_cg_algo(opt_in->ls_cg_algo); // set the cg algorithm using input option

    }


    void RCG::get_search_direction()
    {
        if( iter == 0 || 
            // mani->metric(x1, d1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon())
            
            // FLT_EPSILON      = 1.19209e-07
            // DBL_EPSILON      = 2.22045e-16
            // LDBL_EPSILON     = 1.0842e-19
            mani->metric(x1, d1, gf1) / ngf1 / ngf1 >= -std::sqrt(std::numeric_limits<double>::epsilon() )
            )
        {

            // if(iter > 0 && verbose)
            if(iter > 0 && verbosity == HIGH)
            {
                // d1.brief_print("inside get_search_direction, d1:");
                std::cout << "mani->metric(x1, d1, gf1) / ngf1 / ngf1 = " 
                         << mani->metric(x1, d1, gf1) / ngf1 / ngf1 
                         << std::endl;
            }

            d1 = -1.0 * gf1;

        }
        // else
        // {
        //     d1 = -1.0 * gf1;
        // }
        // d1 = -1.0 * gf1;
    }

    void RCG::update_data()
    {
        // compute_sigma();
        // ManifoldVector Td1, Tgf1;

        // Td1 =mani-> vector_transport(x2, d1);
        double sigmaFR, sigmaPR;

        ManifoldVector Td1, Tgf1;
        ManifoldVector y1;
        Td1 =mani-> vector_transport(x2, d1, x2, d1);  nV ++;
        Tgf1 =mani-> vector_transport(x2, gf1, x2, d1); nV ++;
        // ManifoldVector Td1, Tgf1;
        // ManifoldVector y1;

        if(fabs(mani->metric(x2, gf2, Tgf1)) / mani->metric( x1, gf1, gf1) > 0.1)
        {
            sigma = 0;
        }
        else
        {
            switch(cg_algo)
            {
                case FLETCHER_REEVES:
                    sigma = mani->metric(x2, gf2, gf2) 
                                /mani->metric(x1, gf1, gf1);
                    break;
                case DAI_YUAN:
                    y1 = gf2 - Tgf1;
                    sigma = mani->metric(x1, gf2, gf2) 
                                /mani->metric(x2, Td1, y1);
                    break;
                case POLAK_RIBIERES:
                    sigma = mani->metric(x2, gf2, Tgf1) 
                                /mani->metric(x1, gf1, gf1);
                    break;

                case POLAK_RIBIERES_MOD:
                    y1 = gf2 - Tgf1;
                    sigma = mani->metric(x2, y1, gf2) 
                            /mani->metric(x1, gf1, gf1);
                    sigma = (sigma < 0) ? 0 : sigma;
                    break;

                case FR_PR:
                    sigmaFR = mani->metric(x2, gf2, gf2) 
                                /mani->metric(x1, gf1, gf1);
                    sigmaPR = mani->metric(x2, gf2, Tgf1) 
                                /mani->metric(x1, gf1, gf1);
                    sigma = (sigmaPR < -sigmaFR) ? - sigmaFR : ((sigmaPR > sigmaFR) ? sigmaFR : sigmaPR);
                    break;

                case HAGER_ZHANG:
                    break;

                case HESTENES_STIEFEL:
                    y1 = gf2 - Tgf1;
                    sigma = mani->metric(x1, gf2, y1) 
                                /mani->metric(x2, Td1, y1);
                    break;
                default:
                    std::cout << "Unknown RCG algorithm" << std::endl;
                    sigma = 0;
            }
                    
        }
        // if (cg_algo == FLETCHER_REEVES)
        // {
        //     Td1 =mani-> vector_transport(x2, d1);
        //     sigma = mani->metric(x2, gf2, gf2) 
        //                 /mani->metric(x1, gf1, gf1);
        // }
        // else if (cg_algo == POLAK_RIBIERES ) 
        // {
        //     Td1 =mani-> vector_transport(x2, d1);
        //     sigma = mani->metric(x2, Td1, gf2) 
        //                 /mani->metric(x1, gf1, gf1);
        // }
        // else if (cg_algo == HESTENES_STIEFEL)
        // {
        //     Td1 =mani-> vector_transport(x2, d1);
        //     Tgf1 =mani-> vector_transport(x2, gf1);

        //     y1 = gf2 - Tgf1;
        //     sigma = mani->metric(x1, gf2, y1) 
        //                 /mani->metric(x2, Td1, y1);
        // }

        // else if (cg_algo == DAI_YUAN)
        // {

        //     Td1 =mani-> vector_transport(x2, d1);
        //     Tgf1 =mani-> vector_transport(x2, gf1);
        //     y1 = gf2 - Tgf1;
        //     sigma = mani->metric(x1, gf2, gf2) 
        //                 /mani->metric(x2, Td1, y1);
        // }
        if (sigma == 0)
        {
            d1 = -gf2;
        }
        else 
        {
            d1 = - gf2 + sigma * Td1; // update new direction
        }
        // std::cout << "sigma=" << sigma << std::endl;
    }



    void RCG::set_cg_algo( std::string & name)
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
    
    void RCG::compute_sigma()
    {
        ManifoldVector Td1, Tgf1;
        ManifoldVector y1;
        Td1 =mani-> vector_transport(x2, d1, x2, d1);  nV ++;
        Tgf1 =mani-> vector_transport(x2, gf1, x2, gf1); nV ++;

        if(fabs(mani->metric(x2, gf2, Tgf1)) / mani->metric( x1, gf1, gf1) > 0.1)
        {
            sigma = 0;
        }
        switch(cg_algo)
        {
            case FLETCHER_REEVES:
                sigma = mani->metric(x2, gf2, gf2) 
                            /mani->metric(x1, gf1, gf1);
                break;
            case DAI_YUAN:
                y1 = gf2 - Tgf1;
                sigma = mani->metric(x1, gf2, gf2) 
                            /mani->metric(x2, Td1, y1);
                break;
            case POLAK_RIBIERES:
                sigma = mani->metric(x2, Td1, gf2) 
                        /mani->metric(x1, gf1, gf1);
                break;
            case HAGER_ZHANG:
                break;
            case HESTENES_STIEFEL:
                y1 = gf2 - Tgf1;
                sigma = mani->metric(x1, gf2, y1) 
                            /mani->metric(x2, Td1, y1);
                break;
            default:
                std::cout << "Unknown RCG algorithm" << std::endl;
                sigma = 0;
        }
    }


    // void RCG::update_data()
    // {
    //     ;
    // }
    // void RCG::print_info()
    // {

    // }
} // namespace Module_Optimizer