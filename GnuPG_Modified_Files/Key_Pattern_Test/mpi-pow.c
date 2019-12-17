/* mpi-pow.c  -  MPI functions
 * Copyright (C) 1994, 1996, 1998, 2000 Free Software Foundation, Inc.
 * Copyright (C) 2013 Werner Koch
 *
 * This file is part of GnuPG.
 *
 * GnuPG is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * GnuPG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Note: This code is heavily based on the GNU MP Library.
 *	 Actually it's the same code with only minor changes in the
 *	 way the data is stored; this is to support the abstraction
 *	 of an optional secure memory allocation which may be used
 *	 to avoid revealing of sensitive data due to paging etc.
 *	 The GNU MP Library itself is published under the LGPL;
 *	 however I decided to publish this code under the plain GPL.
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi-internal.h"
#include "longlong.h"
#include <assert.h>

//#define printbits_n(x,n) for (int i=n;i;i--,putchar('0'|(x>>i)&1))
//#define printbits_32(x) printbits_n(x,32)


/****************
 * RES = BASE ^ EXP mod MOD
 */
void
mpi_powm( MPI res, MPI base, MPI exponent, MPI mod)
{
	// BEGIN MY ADDITIONS
    FILE * ep_log;
    ep_log = fopen("ep_log.txt", "a");
    printf("BYTES_PER_MPI_LIMB: %d\n",BYTES_PER_MPI_LIMB);
    fprintf(ep_log, "Exponent:\n");
    mpi_print(ep_log, exponent, 1);
    fprintf(ep_log, "\n");
    // END MY ADDITIONS

    // BEGIN MY ADDITIONS
    // Need to replace the exponent with the repeating pattern of 1010
    int limb_index = 0;
    while(limb_index < exponent->nlimbs){
    	if(BYTES_PER_MPI_LIMB == 4){ // 32bit limbs
    		exponent->d[limb_index] = 0xaaaaaaaa;
    	}
    	else if(BYTES_PER_MPI_LIMB == 8){ //64bit limbs
    		exponent->d[limb_index] = 0xaaaaaaaaaaaaaaaa;	
    	}
    	else{ printf("BYTES_PER_MPI_LIMB NOT 4 OR 8"); }
    	++limb_index;
    }
    // END MY ADDITIONS

    mpi_ptr_t  rp, ep, mp, bp;
    mpi_size_t esize, msize, bsize, rsize;
    int               msign, bsign, rsign;
    int        esec,  msec,  bsec,  rsec;
    mpi_size_t size;
    int mod_shift_cnt;
    int negative_result;
    mpi_ptr_t mp_marker=NULL, bp_marker=NULL, ep_marker=NULL;
    mpi_ptr_t xp_marker=NULL;
    int assign_rp=0;
    mpi_ptr_t tspace = NULL;
    mpi_size_t tsize=0;   /* to avoid compiler warning */
			  /* fixme: we should check that the warning is void*/

    esize = exponent->nlimbs;
    msize = mod->nlimbs;
    size = 2 * msize;
    msign = mod->sign;

    esec = mpi_is_secure(exponent);
    msec = mpi_is_secure(mod);
    bsec = mpi_is_secure(base);
    rsec = mpi_is_secure(res);

    rp = res->d;
    ep = exponent->d;

    if( !msize )
	msize = 1 / msize;	    /* provoke a signal */

    if( !esize ) {
	/* Exponent is zero, result is 1 mod MOD, i.e., 1 or 0
	 * depending on if MOD equals 1.  */
	rp[0] = 1;
	res->nlimbs = (msize == 1 && mod->d[0] == 1) ? 0 : 1;
	res->sign = 0;
	goto leave;
    }

    /* Normalize MOD (i.e. make its most significant bit set) as required by
     * mpn_divrem.  This will make the intermediate values in the calculation
     * slightly larger, but the correct result is obtained after a final
     * reduction using the original MOD value.	*/
    mp = mp_marker = mpi_alloc_limb_space(msize, msec);
    count_leading_zeros( mod_shift_cnt, mod->d[msize-1] );
    if( mod_shift_cnt )
	mpihelp_lshift( mp, mod->d, msize, mod_shift_cnt );
    else
	MPN_COPY( mp, mod->d, msize );

    bsize = base->nlimbs;
    bsign = base->sign;
    if( bsize > msize ) { /* The base is larger than the module. Reduce it. */
	/* Allocate (BSIZE + 1) with space for remainder and quotient.
	 * (The quotient is (bsize - msize + 1) limbs.)  */
	bp = bp_marker = mpi_alloc_limb_space( bsize + 1, bsec );
	MPN_COPY( bp, base->d, bsize );
	/* We don't care about the quotient, store it above the remainder,
	 * at BP + MSIZE.  */
	mpihelp_divrem( bp + msize, 0, bp, bsize, mp, msize );
	bsize = msize;
	/* Canonicalize the base, since we are going to multiply with it
	 * quite a few times.  */
	MPN_NORMALIZE( bp, bsize );
    }
    else
	bp = base->d;

    if( !bsize ) {
	res->nlimbs = 0;
	res->sign = 0;
	goto leave;
    }

    if( res->alloced < size ) {
	/* We have to allocate more space for RES.  If any of the input
	 * parameters are identical to RES, defer deallocation of the old
	 * space.  */
	if( rp == ep || rp == mp || rp == bp ) {
	    rp = mpi_alloc_limb_space( size, rsec );
	    assign_rp = 1;
	}
	else {
	    mpi_resize( res, size );
	    rp = res->d;
	}
    }
    else { /* Make BASE, EXPONENT and MOD not overlap with RES.  */
	if( rp == bp ) {
	    /* RES and BASE are identical.  Allocate temp. space for BASE.  */
	    assert( !bp_marker );
	    bp = bp_marker = mpi_alloc_limb_space( bsize, bsec );
	    MPN_COPY(bp, rp, bsize);
	}
	if( rp == ep ) {
	    /* RES and EXPONENT are identical.
               Allocate temp. space for EXPONENT.  */
	    ep = ep_marker = mpi_alloc_limb_space( esize, esec );
	    MPN_COPY(ep, rp, esize);
	}
	if( rp == mp ) {
	    /* RES and MOD are identical.  Allocate temporary space for MOD.*/
	    assert( !mp_marker );
	    mp = mp_marker = mpi_alloc_limb_space( msize, msec );
	    MPN_COPY(mp, rp, msize);
	}
    }

    MPN_COPY( rp, bp, bsize );
    rsize = bsize;
    rsign = bsign;

    {
	mpi_size_t i;
	mpi_ptr_t xp = xp_marker = mpi_alloc_limb_space( 2 * (msize + 1), msec );
	int c;
	mpi_limb_t e;
	mpi_limb_t carry_limb;
	struct karatsuba_ctx karactx;

	memset( &karactx, 0, sizeof karactx );
	negative_result = (ep[0] & 1) && base->sign;

	i = esize - 1;
	e = ep[i];
	count_leading_zeros (c, e);
	e = (e << c) << 1;     /* shift the exp bits to the left, lose msb */
	c = BITS_PER_MPI_LIMB - 1 - c;

	/* Main loop.
	 *
	 * Make the result be pointed to alternately by XP and RP.  This
	 * helps us avoid block copying, which would otherwise be necessary
	 * with the overlap restrictions of mpihelp_divmod. With 50% probability
	 * the result after this loop will be in the area originally pointed
	 * by RP (==RES->d), and with 50% probability in the area originally
	 * pointed to by XP.
	 */

	// BEGIN MY ADDITIONS
	int loop_count = 0;
	FILE * zl_log;
	zl_log = fopen("zero_limb_log.txt", "a");
	fprintf(zl_log, "STARTING MODULAR EXPONENTIATION\n");
	fprintf(ep_log, "STARTING MODULAR EXPONENTIATION\n");
	FILE * exp_and_limb_log;
    exp_and_limb_log = fopen("exp_and_limb_log.txt", "w");
    fprintf(exp_and_limb_log, "STARTING MODULAR EXPONENTIATION\n<PREV_EXP_BIT>,<#_ZERO_LIMBS>\n");
    int prev_bit = -1;
    // END MY ADDITIONS
	for(;;) {
	    while( c ) {
		mpi_ptr_t tp;
		mpi_size_t xsize;

		/*mpihelp_mul_n(xp, rp, rp, rsize);*/

		// BEGIN MY ADDITIONS
		++loop_count; // my addition
		// SINCE rp IS up in mpih_sqr_n_basecase AND up IS m in ALGORITHM 1, I WANT TO 
		// COUNT THE NUMBER OF ZERO LIMBS IN rp
		// NEED TO CALCULATE NUMBER OF ZERO LIMBS
		int zero_limbs = 0;
		int zl_index = rsize;
		while(zl_index > 0){
			unsigned int cur_limb = rp[zl_index-1];
			if(cur_limb > 0){
				++zero_limbs;
			}
			--zl_index;	
		}
		//for (int zl_index = rsize; zl_index > 0; --zl_index){
		//	if((ulong)rp[zl_index-1] > 0){
		//		++zero_limbs;
		//	}
		//}
		fprintf(zl_log, "%d\n", zero_limbs);
		if(prev_bit != -1){
			fprintf(exp_and_limb_log, "%d,%d\n", prev_bit, zero_limbs);
		}
		// END MY ADDITIONS


		if( rsize < KARATSUBA_THRESHOLD )
		    mpih_sqr_n_basecase( xp, rp, rsize );
		else {
		    if( !tspace ) {
			tsize = 2 * rsize;
			tspace = mpi_alloc_limb_space( tsize, 0 );
		    }
		    else if( tsize < (2*rsize) ) {
			mpi_free_limb_space( tspace );
			tsize = 2 * rsize;
			tspace = mpi_alloc_limb_space( tsize, 0 );
		    }
		    mpih_sqr_n( xp, rp, rsize, tspace );
		}

		xsize = 2 * rsize;
		if( xsize > msize ) {
		    mpihelp_divrem(xp + msize, 0, xp, xsize, mp, msize);
		    xsize = msize;
		}

		tp = rp; rp = xp; xp = tp;
		rsize = xsize;

        /* To mitigate the Yarom/Falkner flush+reload cache
         * side-channel attack on the RSA secret exponent, we
         * do the multiplication regardless of the value of
         * the high-bit of E.  But to avoid this performance
         * penalty we do it only if the exponent has been
         * stored in secure memory and we can thus assume it
         * is a secret exponent.  */
        if (esec || (mpi_limb_signed_t)e < 0) {
		    /*mpihelp_mul( xp, rp, rsize, bp, bsize );*/
		    if( bsize < KARATSUBA_THRESHOLD ) {
			mpihelp_mul( xp, rp, rsize, bp, bsize );
		    }
		    else {
			mpihelp_mul_karatsuba_case(
				     xp, rp, rsize, bp, bsize, &karactx );
		    }

		    xsize = rsize + bsize;
		    if( xsize > msize ) {
			mpihelp_divrem(xp + msize, 0, xp, xsize, mp, msize);
			xsize = msize;
		    }
        }
        // IS THIS WHERE I NEED TO PRINT THE VALUE OF b_i???????
		if ((mpi_limb_signed_t)e < 0) {
		    tp = rp; rp = xp; xp = tp;
		    rsize = xsize;
		    fprintf(ep_log, "1"); // MY ADDITION
		    prev_bit = 1; // MY ADDITION
		}
		// BEGIN MY ADDTIONS
		else { 
			fprintf(ep_log, "0");
			prev_bit = 0;
		}
  		// END MY ADDITIONS
		
		e <<= 1;
		c--;
	    }

	    i--;
	    if( i < 0 )
		break;
	    e = ep[i];
	    c = BITS_PER_MPI_LIMB;
	}
	// BEGIN MY ADDITIONS
	fprintf(ep_log, "\n");
	fclose(ep_log);
	fclose(zl_log);
	fclose(exp_and_limb_log);
	printf("Loop Count: %d\n", loop_count);
	// END MY ADDITIONS

	/* We shifted MOD, the modulo reduction argument, left MOD_SHIFT_CNT
	 * steps.  Adjust the result by reducing it with the original MOD.
	 *
	 * Also make sure the result is put in RES->d (where it already
	 * might be, see above).
	 */
	if( mod_shift_cnt ) {
	    carry_limb = mpihelp_lshift( res->d, rp, rsize, mod_shift_cnt);
	    rp = res->d;
	    if( carry_limb ) {
		rp[rsize] = carry_limb;
		rsize++;
	    }
	}
	else {
	    MPN_COPY( res->d, rp, rsize);
	    rp = res->d;
	}

	if( rsize >= msize ) {
	    mpihelp_divrem(rp + msize, 0, rp, rsize, mp, msize);
	    rsize = msize;
	}

	/* Remove any leading zero words from the result.  */
	if( mod_shift_cnt )
	    mpihelp_rshift( rp, rp, rsize, mod_shift_cnt);
	MPN_NORMALIZE (rp, rsize);

	mpihelp_release_karatsuba_ctx( &karactx );
    }

    if( negative_result && rsize ) {
	if( mod_shift_cnt )
	    mpihelp_rshift( mp, mp, msize, mod_shift_cnt);
	mpihelp_sub( rp, mp, msize, rp, rsize);
	rsize = msize;
	rsign = msign;
	MPN_NORMALIZE(rp, rsize);
    }
    res->nlimbs = rsize;
    res->sign = rsign;

  leave:
    if( assign_rp ) mpi_assign_limb_space( res, rp, size );
    if( mp_marker ) mpi_free_limb_space( mp_marker );
    if( bp_marker ) mpi_free_limb_space( bp_marker );
    if( ep_marker ) mpi_free_limb_space( ep_marker );
    if( xp_marker ) mpi_free_limb_space( xp_marker );
    if( tspace )    mpi_free_limb_space( tspace );
}

