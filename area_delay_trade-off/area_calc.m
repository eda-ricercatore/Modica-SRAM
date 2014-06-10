%{
This is written by Zhiyang Ong (zhiyango@usc.edu;
6004 9194 12) and Andrew Mattheisen (mattheis@usc.edu;
2134 5147 11) for EE 577B, SRAM Project - Part 2, Questions
1 and 2

Area calculation for the SRAM
@param num_rows is the number of rows in the SRAM block
@param num_cols is the number of columns in the SRAM block
@param cur_area is the current area in the SRAM block
@return the total area for the SRAM design


The area for the chip is calculated by the square of its
longer dimension.

Thus, if the SRAM block is rectangular, and the square of
its longer dimension can fit all of the remaining SRAM
circuitry, this square shall be the total chip area of the
SRAM design.

Else, do the following.
Add the remaining circuitry to the SRAM design on the side
with the shorter dimension.
The circuit and the SRAM block are abutted on their longer
edges/dimensions, so that the longer dimension of the SRAM
block will be that for the SRAM design.
This is because the vertical and horizontal dimensions of
the SRAM block are much greater than that of the remaining
circuitry of the SRAM.


Let the dimensions of the SRAM block be x and y, where x<y

Let the dimensions of the remaining circuitry in the SRAM
design be a and b, where a<b

Since a<x and b<y, the area for the SRAM design will be
(x+a) by (y)


Doing otherwise will definitely increase the area of the
SRAM design, since the area of the SRAM block is y^2 and
the area of the SRAM would be (y+a)^2 or (y+b)^2.

All the considered SRAM blocks are rectangular since the
SRAM must be word-addressable.

Thus, either dimensions of the SRAM block cannot be an odd
multiple of bytes, or halfwords
%}

function [t_area] = area_calc(num_rows, num_cols, cur_area)
%{
The lengths for num_rows, num_cols, and cur_area should
be the same, since n designs will be compared; n refers
to the length of these arrays.
%}

    if length(num_rows) ~=  length(num_cols)
        error('Length of num_rows and num_cols do not match')
    end
    
    if length(num_rows) ~=  length(cur_area)
        error('Length of num_rows and cur_area do not match')
    end


    for i=1:length(num_rows)
        % Longer dimension of the SRAM block
        longer_dim(i)=max(num_rows(i), num_cols(i));
        % Shorter dimension of the SRAM block
        shorter_dim(i)=min(num_rows(i), num_cols(i));
        
        %{
            Total area of SRAM block
            = (length of longer dimension)^2 = y^2
        %}
        longer_dim_sq(i)= longer_dim(i)^2;
        %{
            Actual area of the SRAM block
            = length of the block * height of the block=x*y
        %}
        actual_area(i)=num_rows(i)*num_cols(i);
        % Determine the leftover/remaining space: y*(y-x)
        remainder(i)=longer_dim(i)*(longer_dim(i)-shorter_dim(i));
        %{
            Can the remaining circuit fit into the empty
            space?
            
            This is done because if the longer dimension
            was halved, it would not work for the 1 * 2^{15}
        
            If not, abut the remaining circuit to the SRAM
            block
        %}
        if cur_area(i) < remainder(i)
            %{
                The remaining circuit can fit into the
                leftover space
            %}
            t_area(i) = longer_dim_sq(i);
        else
            %{
                Abut the remaining circuit to the SRAM 
                block
                
                Since the exact dimensions of the remaining
                circuit in the SRAM design can only be known
                after elaborate estimation of its delay to
                size its gates, the authors assume that the
                remaining circuit can be evenly spread
                around the edges of the SRAM block to obtain
                the final dimensions of the SRAM design
            
                This assumption allows us to sum the area of
                the SRAM block and the remaining circuit
                to obtain the final dimensions of the SRAM
                design
            %}
            t_area(i) = longer_dim_sq(i)+cur_area(i);
        end
        
%{
Note that I should have explicity considered the area taken 
up by wiring in the SRAM design.

However, in the remaining circuit, most of the area is taken
up by the larger of the row or column decoder.

Since the height of the row/column decoder is given by twice
the wiring pitch needed to route to all rows/columns, and 
its length is much less than that of the SRAM block, wiring
does not have a prominent presence in our SRAM designs.

This is because 40 \lambda is greater than twice the wiring
pitch per row/column, and the number of stages required in
the remaining circuit do not exceed 40 \lambda * the number
of columns/rows.
%}
    end



