ONE   = 1.0e0
TWO   = 2.0e0
THREE = 3.0e0
HALF  = 0.5e0
Sqrt3 = Math.sqrt(3.0e0)
Sqrt2 = Math.sqrt(2.0e0)

#EMPTY =  1 - 1
#UP    =  2 - 1
#DOWN  =  3 - 1
#FULL  =  4 - 1
#UNITY =  1 - 1
#REPR_STATE = [ 1, 2, 2, 1 ]
#A     =  2 - 1
#C     =  3 - 1
#AAS   =  4 - 1
#CAS   =  5 - 1
#CAT   =  6 - 1
#CCS   =  7 - 1
#CAA   =  8 - 1
#CCA   =  9 - 1
#CCAA  = 10 - 1
#REPR_OP = [ 1, 2, 2, 1, 1, 3, 1, 2, 2, 1 ]

EMPTY  =  1
UP     =  2
DOWN   =  3
FULL   =  4

UNITY  =  1
A      =  2
C      =  3
AAS    =  4
CAS    =  5
CAT    =  6
CCS    =  7
CAA    =  8
CCA    =  9
CCAA   = 10
# addtional begin
ACS    = 11
ACT    = 12
CASCAS = 13
# addtional end
REPR_STATE = [ 0, 1, 2, 2, 1 ]
REPR_OP = [ 0, 1, 2, 2, 1, 1, 3, 1, 2, 2, 1, 1, 3, 1 ]

TBL_WFSTR = [ "NULL", "EMPTY", "UP", "DOWN", "FULL" ]
TBL_OPSTR = [ "NULL", "UNITY", "A", "C", "AAS", "CAS", "CAT", "CCS", "CAA", "CCA", "CCAA", "ACS", "ACT", "CASCAS" ]

def query_stepval( step )
    if step == 0 then
        return EMPTY
    elsif step == 1 then
        return UP
    elsif step == 2 then
        return DOWN
    elsif step == 3 then
        return FULL
    else
        $stderr.printf( "Error: query_stepval: 0<= argument < 4 : %s\n", step )
        exit
    end
end
