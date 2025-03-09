        program test_va_arg_fortran_2_c
        implicit none
        character  sp*80
        integer test_va_arg_f3, val_i,i

        sp="hello world"//char(0)
        val_i = 666
        i=test_va_arg_f3(val_i,sp)
        end
