        integer m_pix, m_mdls
        parameter( m_pix = 131072 )
        parameter( m_mdls = 8 )

        integer n_comp, m_crrnt, p_model, model_( m_mdls )
        real m_e( m_mdls ), m_rat( m_mdls )
        real norm( m_mdls ), param( 3, m_mdls ), nh( m_mdls )
        real z1, nh_g
        common / no1_cup / n_comp, m_crrnt, p_model, model_, m_e, m_rat,
     &                                         norm, param, nh, z1, nh_g

        real e_in( m_pix ), f_in( m_pix ), e_old( m_mdls )
        integer beg_m( m_mdls ), fin_m( m_mdls ), j_old( m_mdls )
        integer highest
        logical outob( m_mdls )
        common / no2_cup / e_in, f_in, e_old, beg_m, fin_m,
     &                                             j_old, highest, outob

        character*128 mf_name( m_mdls )
        common / no3_cup / mf_name

        real bg_rate
        common / no7_cup / bg_rate
