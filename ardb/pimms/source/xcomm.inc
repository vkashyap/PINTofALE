	Integer m_Par, l_Buffer
	Integer m_Comm, l_Comm, m_Buffer, n_Res
	Parameter( m_Par = 32 )
	Parameter( l_Buffer = 256 )
	Parameter( m_Comm = 128, l_Comm = 16 )
        parameter( m_Buffer = 64 )
	Parameter( n_Res = 2 )

        Character*( l_buffer ) current_Line
        common / fci_00c / current_Line
        character*( l_comm ) list( m_comm ), prompt
        common / fci_01c / list, prompt
        integer n_comm, order( m_comm ), length
        common / fci_01i / n_comm, order, length
        logical empty, from_file
        common / fci_01l / empty, from_file
        character*( l_buffer ) buff_line( m_buffer )
        common / fci_02c / buff_line
        integer n_buffer, file_unit
        common / fci_02i / n_buffer, file_unit
        integer fci_reserved( n_res )
        common / fci_03i / fci_reserved
        character*( l_buffer ) p_char
        common / fci_04c / p_char
        integer pc_b( m_par ), pc_f( m_par ), np_c
        common / fci_04i / pc_b, pc_f, np_c
        real p_float( m_par )
        integer p_int( m_par ), np_f, np_i, n_list( m_par * 2 ), np_n
        common / fci_05 / p_float, p_int, np_f, np_i, n_list, np_n
        integer np_total, p_kind( m_par * 3 ), po_c( m_par )
        integer po_f( m_par ), po_i( m_par ), po_n( m_par )
        common / fci_05b / np_total, p_kind, po_c, po_f, po_i, po_n
