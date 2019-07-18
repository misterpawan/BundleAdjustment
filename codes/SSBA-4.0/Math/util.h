// This file contains the utility functions used in the code

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80
int file_row_count ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  string line;
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

void cc_header_read ( string prefix, int &ncc, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_HEADER_READ reads header information about a sparse matrix in CC format.
//
//  Discussion:
//
//    Three files are presumed to exist:
//    * prefix_icc.txt contains NCC ICC values;
//    * prefix_ccc.txt contains N+1 CCC values;
//    * prefix_acc.txt contains NCC ACC values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, a common prefix for the filenames.
//
//    Output, int &NCC, the number of CC elements.
//
//    Output, int &N, the number of columns in the matrix.
//
{
  string filename_ccc;
  string filename_icc;

  filename_icc = prefix + "_row.txt";
  ncc = file_row_count ( filename_icc );

  filename_ccc = prefix + "_col.txt";
  n = file_row_count ( filename_ccc ) - 1;

  return;
}

void r8vec_data_read ( string input_filename, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DATA_READ reads the data from an R8VEC file.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double TABLE[N], the data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  int lchar;
  string line;
  double x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8VEC_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    table[j] = atof ( line.c_str ( ) );
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

//****************************************************************************80

void i4vec_data_read ( string input_filename, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DATA_READ reads data from an I4VEC file.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly one value.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int TABLE[N], the data.
//
{
  ifstream input;
  int i;
  int j;
  int l;
  string line;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4VEC_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    table[j] = atoi ( line.c_str ( ) );
    j = j + 1;
  }

  input.close ( );

  return;
}

//****************************************************************************80

void cc_data_read ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_DATA_READ reads data about a sparse matrix in CC format.
//
//  Discussion:
//
//    Three files are presumed to exist:
//    * prefix_icc.txt contains NCC ICC values;
//    * prefix_ccc.txt contains N+1 CCC values;
//    * prefix_acc.txt contains NCC ACC values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, a common prefix for the filenames.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int N, the number of columns in the matrix.
//
//    Output, int ICC[NCC], the CC rows.
//
//    Output, int CCC[N+1], the compressed CC columns.
//
//    Output, double ACC[NCC], the CC values.
//
{
  string filename_acc;
  string filename_ccc;
  string filename_icc;

  filename_icc = prefix + "_row.txt";
  i4vec_data_read ( filename_icc, ncc, icc );

  filename_ccc = prefix + "_col.txt";
  i4vec_data_read ( filename_ccc, n + 1, ccc );

  filename_acc = prefix + "_val.txt";
  r8vec_data_read ( filename_acc, ncc, acc );

  return;
}
