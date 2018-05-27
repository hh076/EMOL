
module EmolUtil

  def transpose(mat, m = nil)
    m ||= mat[0].size
    (0..m-1).collect { |i| mat.collect { |a| a[i] } }
  end

  def print_2dim_ary(ary, var, print_row_size, fmt)
    row_size = ary[0].size; clm_size = ary.size

    printf("\n  %8s \n", var)
    if row_size <= print_row_size
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], 0, row_size-2, " \n", ary[i][row_size-1], fmt)
      end
    else
      from = 0
      while row_size - from > print_row_size
        printf("\n")
        for i in 0..clm_size-1
          printf("%5d   ", i+1)
          print_line(ary[i], from, print_row_size-2, " \n", ary[i][from+print_row_size-1], fmt)
        end
        from += print_row_size
      end

      printf("\n")
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], from, row_size-2-from, " \n", ary[i][row_size-1], fmt)
      end
    end
  end

  def print_2dim_ary_with_value(ary, value, var, print_row_size, fmt)
    row_size = ary[0].size; clm_size = ary.size
#   row_size = value.size; clm_size = ary.size

    printf("\n  %20s \n", var)
    if row_size <= print_row_size
      print("        ")
      print_line(value, 0, row_size-2, " \n", value[row_size-1], fmt)
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], 0, row_size-2, " \n", ary[i][row_size-1], fmt)
      end
    else
      from = 0
      while row_size - from > print_row_size
        printf("\n        ")
        print_line(value, from, print_row_size-2, " \n", value[from+print_row_size-1], fmt)
        for i in 0..clm_size-1
          printf("%5d   ", i+1)
          print_line(ary[i], from, print_row_size-2, " \n", ary[i][from+print_row_size-1], fmt)
        end
        from += print_row_size
      end

      printf("\n        ")
      print_line(value, from, row_size-2-from, " \n", value[row_size-1], fmt)
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], from, row_size-2-from, " \n", ary[i][row_size-1], fmt)
      end
    end
  end

  def print_2dim_ary_with_value_moname(ary, value, bs_name, var, print_row_size, fmt)
    row_size = ary[0].size; clm_size = ary.size
#   row_size = value.size; clm_size = ary.size

    printf("\n  %20s \n", var)
    if row_size <= print_row_size
      printf("               ")
      ( 1..row_size ).each do | i |
          printf( "         %5d", i )
      end
      printf("\n                 ")

      print_line(value, 0, row_size-2, " \n", value[row_size-1], fmt)
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        printf("%6s   ", bs_name[i])
        print_line(ary[i], 0, row_size-2, " \n", ary[i][row_size-1], fmt)
      end
    else
      from = 0
#     from = 1
      while row_size - from > print_row_size
        printf("               ")
      # print_line(num, from, print_row_size-2, " \n", value[from+print_row_size-1], "       %5d")
        ( from..(from + print_row_size - 1) ).each do | i |
            printf( "         %5d", i + 1 )
        end
        printf("\n")
        printf("                 ")
        print_line(value, from, print_row_size-2, " \n", value[from+print_row_size-1], fmt)
        for i in 0..clm_size-1
          printf("%5d   ", i+1)
          printf("%6s   ", bs_name[i])
          print_line(ary[i], from, print_row_size-2, " \n", ary[i][from+print_row_size-1], fmt)
        end
        from += print_row_size
        printf("\n")
      end

      printf("               ")
      ( from..row_size-1 ).each do | i |
          printf( "         %5d", i + 1 )
      end
 #    print_line(num, from, print_row_size-2-from, " \n", value[from+print_row_size-1], "      %5d")
      printf("\n                 ")
      print_line(value, from, row_size-2-from, " \n", value[row_size-1], fmt)
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        printf("%6s   ", bs_name[i])
        print_line(ary[i], from, row_size-2-from, " \n", ary[i][row_size-1], fmt)
      end
    end
  end

  def print_2dim_ary_tri(ary, var, print_row_size, fmt)
    row_size = ary[0].size; clm_size = ary.size

    printf("\n  %8s \n", var)
    if row_size <= print_row_size
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], 0, row_size-2, " \n", ary[i][row_size-1], fmt)
      end
    else
      from = 0
      while row_size - from > print_row_size
        printf("\n")
        for i in 0..clm_size-1
          printf("%5d   ", i+1)
          print_line(ary[i], from, print_row_size-2, " \n", ary[i][from+print_row_size-1], fmt)
        end
        from += print_row_size
      end

      printf("\n")
      for i in 0..clm_size-1
        printf("%5d   ", i+1)
        print_line(ary[i], from, row_size-2-from, " \n", ary[i][row_size-1], fmt)
      end
    end
  end

  def print_ary(ary, var, print_row_size, fmt)
    printf("\n       %8s = [ ", var)
    if ary.size <= print_row_size
      print_line(ary, 0, ary.size-2, " ] \n", ary[ary.size-1], fmt)
    else
      from = 0
      while ary.size - from > print_row_size
        printf("\n            ")
        print_line(ary, from, print_row_size-2, ", \\", ary[from+print_row_size-1], fmt)
        from += print_row_size
      end
      printf("\n            ")
      print_line(ary, from, ary.size-2-from, " ] \n", ary[ary.size-1], fmt)
    end
  end

  def print_ary_tri(ary, var, print_row_size, fmt)
    printf("\n       %8s = [ \n", var)
    row = 1; from = 0
    while ary.size - from >= row
      printf("%5d       ", row)
      if row <= print_row_size
        print_line(ary, from, row-2, " \n", ary[from+row-1], fmt)
        from += row; row += 1
      else
        last = from + row -1
        while last - from + 1 > print_row_size
          print_line(ary, from, print_row_size-2, ", \\ \n", ary[from+print_row_size-1], fmt)
          printf("            ")
          from += print_row_size
        end
        print_line(ary, from, last-from-1, " \n", ary[last], fmt)
        from = last + 1; row += 1
      end
    end
    printf("             ] \n")
  end

  def print_line(a, from, to, br, a_last, fmt)
    for i in 0..to do
#     printf(fmt + ", ", a[from + i])    
      printf(fmt + "  ", a[from + i])    
    end
    printf(fmt + br, a_last)     
  end

  module_function :transpose
  module_function :print_ary
  module_function :print_ary_tri
  module_function :print_2dim_ary
  module_function :print_2dim_ary_with_value
  module_function :print_2dim_ary_with_value_moname
  module_function :print_2dim_ary_tri
  module_function :print_line
end
