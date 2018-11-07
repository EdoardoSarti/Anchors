/////////////////////////////////////////////////////////////////////////
//!
//  This file is part of the AlignMe program
//
//  Copyright (C) 2010 by Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//  AlignMe@rzg.mpg.de
//
//  AlignMe is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  AlignMe is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//!
//!
//!  anchor.h
//!
//! Anchors allow to define residues that are kept aligned.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////



#ifndef ANCHOR_H_
#define ANCHOR_H_

#include "matrix.t.h"
#include "dynamic_programing_matrix_element.h"

void
BuildAnchorMatrix( std::istream &STREAM, Matrix< DynamicProgrammingMatrixElement> &MATRIX, const double &PENALTY = -1.0)
{
	// anchors need to be given in usual vector notation
	double
		shift_penalty;
	size_t
		nr, ii, jj;

	STREAM >> nr;
	for( size_t r = 0; r < nr; ++r)
	{
		STREAM >> ii >> jj >> shift_penalty;
		std::cout << "anchor: " << ii << " " << jj << " " << shift_penalty << std::endl;

//		ii -= 1;
//		jj -= 1;
		ii += 1;  // first residue id is 0
		jj += 1;
		shift_penalty *= -1;  // using positive values in file!

		if( ii >= MATRIX.GetNumberOfRows() || jj >= MATRIX.GetNumberOfColumns())
		{
			std::cout << "provide pairs of values that match the lengths of the sequences:" << "\n";
			std::cout << ii << " should be smaller: " << MATRIX.GetNumberOfRows() << "\n";
			std::cout << jj << " should be smaller: " << MATRIX.GetNumberOfColumns() << "\nbye\n\n";
			exit( -1);
		}
		//std::cout << "anchors: " << ii << "  " << jj << "  " << shift_penalty << "\n";

		for( size_t i = 0; i < MATRIX.GetNumberOfRows(); ++i)
		{
			for( size_t j = 0; j < MATRIX.GetNumberOfColumns(); ++j)
			{
//				if( i < ii && j == jj)
				if( i == ii && j < jj)
				{
					for( size_t k = 0; k < 6; k++)
					{
						MATRIX( i, j).AddValue( k, shift_penalty);
					}
				}
//				else if( i == ii && j < jj)
				else if( i < ii && j == jj)
				{
					for( size_t k = 0; k < 3; k++)
					{
						MATRIX( i, j).AddValue( k, shift_penalty);
					}
					for( size_t k = 6; k < 9; k++)
					{
						MATRIX( i, j).AddValue( k, shift_penalty);
					}
				}
//				else if( i == ii && j > jj)
				else if( i > ii && j == jj)
				{
					for( size_t k = 0; k < 3; k++)
                                        {
                                                MATRIX( i, j).AddValue( k, -shift_penalty);
                                        }
					for( size_t k = 6; k < 9; k++)
                                        {
                                                MATRIX( i, j).AddValue( k, -shift_penalty);
                                        }
				}
//				else if( i > ii && j == jj)
				else if( i == ii && j > jj)
                                {
                                        for( size_t k = 0; k < 6; k++)
                                        {
                                                MATRIX( i, j).AddValue( k, -shift_penalty);
                                        }
				}
				else if( i == ii && j == jj)
				{
					for( size_t k = 3; k < 9; k++)
					{
						MATRIX( i, j).AddValue( k, -shift_penalty);
					}
                                }
			}
		}
	}
}



#endif /* ANCHOR_H_ */
