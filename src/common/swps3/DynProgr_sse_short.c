/** \file DynProgr_sse_short.c
 *
 * Profile generation and alignment for packed short vectors on SSE2.
 */
/*
 * Copyright (c) 2007-2008 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "DynProgr_sse_short.h"
#include "Page_size.h"
#include <stdio.h>
#include <float.h>

#define PAGE_ALIGN(x) (((size_t)(x)+getPageSize()-1)&~(getPageSize()-1))
EXPORT ProfileShort * swps3_createProfileShortSSE(const char * query,
		int queryLen, SMatrix matrix) {
	/* int segLen  = ((queryLen+7)/8 + 1) & ~1; */
	int segLen = (queryLen + 7) / 8;
	int i, j, k;
	uint16_t bias = 0;
	uint16_t * pprofile;
	ProfileShort * profile = malloc(
			sizeof(ProfileShort)
					+ ((segLen * MATRIX_DIM + 1) & ~(0x1)) * sizeof(__m128i )
					+ segLen * 3 * sizeof(__m128i ) + 64 + 2 * getPageSize());

	profile->loadOpt = (__m128i *) ((size_t)(profile->data + 15) & ~(0xf));
	profile->storeOpt = profile->loadOpt + segLen;
	profile->rD = profile->storeOpt + segLen;
	profile->profile = (__m128i *) PAGE_ALIGN(profile->rD + segLen);

	/* Init the profile */
	profile->len = queryLen;
	/* Init the byte profile */
	for (i = 0; i < MATRIX_DIM; i++)
		for (j = 0; j < MATRIX_DIM; j++)
			if (bias < -matrix[i * MATRIX_DIM + j])
				bias = -matrix[i * MATRIX_DIM + j];
	pprofile = (uint16_t*) profile->profile;
	/* TODO remove bias? */
	bias = 0;

	for (i = 0; i < MATRIX_DIM; i++) {
		for (j = 0; j < segLen; j++) {
			for (k = 0; k < 8; k++) {
				if (j + k * segLen < queryLen) {
					char queryChar = query[j + k * segLen];
					*(pprofile++) = matrix[queryChar * MATRIX_DIM + i] + bias;
				} else {
					*(pprofile++) = 0;
				}
			}
		}
	}
	profile->bias = bias;
	return profile;
}

EXPORT double swps3_alignmentShortSSE_lin(ProfileShort * query, const char * db,
		int dbLen, Options * options) {

	/**********************************************************************
	 * This version of the code implements the idea presented in
	 *
	 ***********************************************************************
	 * Striped Smith-Waterman speeds database searches six times over other
	 * SIMD implementations
	 *
	 * Michael Farrar, Bioinformatics, 23(2), pp. 156-161, 2007
	 **********************************************************************/

	int i, j;

	uint16_t MaxScore = 0x8000;
	int segLength = (query->len + 7) / 8; /* the segment length */

	__m128i * loadOpt = query->loadOpt;
	__m128i * storeOpt = query->storeOpt;
	__m128i * current_profile;
	__m128i * swap;

	__m128i vMinimums = _mm_set1_epi16(0x8000);

	__m128i vDelFixed = _mm_set1_epi16(options->gapOpen);
	__m128i vBias = _mm_set1_epi16(query->bias);

	__m128i vMaxScore = vMinimums; /* vMaxScore = [0,0] */

	__m128i vProfile = vMinimums; /* the score profile */
	__m128i vStoreOpt; /* the new optimal score */
	__m128i vRD; /* the new row deletion score */
	__m128i vCD; /* the column deletion score */
	__m128i vTmp;

#ifdef PY_DEBUG
	printf("Linear alignment!\n");
#endif

	/* initialize the other arrays used for the dynProg code */
	/*********************************************************/
	for (i = 0; LIKELY(i < segLength); i++) {
		_mm_store_si128(loadOpt + i, vMinimums);
		_mm_store_si128(storeOpt + i, vMinimums);
	}

	/* looping through all the columns */
	/***********************************/
	for (j = 0; LIKELY(j < dbLen); j++) {

		/* compute the opt and cd score depending on the previous column */
		/*******************************************************************/
		/* set the column deletion score to zero, has to be fixed later on */
		vCD = vMinimums;

		/* set the opt score to the elements computed in the previous column */
		/* set the low of storeOpt to MaxS[j] */
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
		vStoreOpt = _mm_slli_si128(vStoreOpt, 2);
		vStoreOpt = _mm_insert_epi16(vStoreOpt, (int) 0x8000, 0);

		/* compute the current profile, depending on the character in s2 */
		/*****************************************************************/

		current_profile = query->profile + db[j] * segLength;
		/* swap the old optimal score with the new one */
		/***********************************************/
		swap = storeOpt;
		storeOpt = loadOpt;
		loadOpt = swap;

		/* main loop computing the max, precomputing etc. */
		/**************************************************/
		for (i = 0; LIKELY(i < segLength); i++) {

			vTmp = _mm_load_si128(loadOpt + i);
			vRD = _mm_adds_epi16(vTmp, vDelFixed);

			/* load the current profile */
			/*vProfile = _mm_movpi64_epi64(current_profile[i]);*/
			/*vProfile = _mm_loadl_epi64((__m128i*)(current_profile+i));*/
			/*#if (defined _WIN32 || defined __WIN32__)*/
			vProfile = _mm_load_si128(current_profile + i);
			/*#else
			 __asm__("MOVDQA (%1),%0" : "=x" (vProfile) : "r" (current_profile+i));
			 #endif*/
			/*vProfile = _mm_unpacklo_epi8(vProfile, _mm_xor_si128(vProfile,vProfile));*/
			vProfile = _mm_subs_epi16(vProfile, vBias);

			/* add the profile the prev. opt */
			vStoreOpt = _mm_adds_epi16(vStoreOpt, vProfile);

			/* update the maxscore found so far */
			vMaxScore = _mm_max_epi16(vMaxScore, vStoreOpt);

			/* compute the correct opt score of the cell */
			vStoreOpt = _mm_max_epi16(vStoreOpt, vCD);
			vStoreOpt = _mm_max_epi16(vStoreOpt, vRD);

			/* store the opt score of the cell */
			_mm_store_si128(storeOpt + i, vStoreOpt);

			/* precompute cd for next iteration */
			vCD = _mm_adds_epi16(vStoreOpt, vDelFixed);

			/* load precomputed opt for next iteration */
			vStoreOpt = vTmp;
		}

		for (i = 0; LIKELY(i < 8); ++i) {
			int k;
			/* compute the gap extend penalty for the current cell */
			vCD = _mm_slli_si128(vCD, 2);
			vCD = _mm_insert_epi16(vCD, 0x8000, 0);

			for (k = 0; LIKELY(k < segLength); ++k) {
				/* compute the current optimal value of the cell */
				vTmp = _mm_load_si128(storeOpt + k);
				vStoreOpt = _mm_max_epi16(vTmp, vCD);
				_mm_store_si128(storeOpt + k, vStoreOpt);

				if (UNLIKELY(
						!_mm_movemask_epi8(_mm_cmpgt_epi16(vStoreOpt, vTmp))))
					goto shortcut;

				/* precompute the scores for the next cell */
				vCD = _mm_adds_epi16(vStoreOpt, vDelFixed);
			}
		}
		shortcut:
		/* store the new MaxScore for the next line block */
		/**************************************************/

		/* store the element of storeOpt in MaxS */
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
	}
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 8));
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 4));
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 2));
	MaxScore = _mm_extract_epi16(vMaxScore, 0);
	if (MaxScore == 0x7fff) {
		return DBL_MAX;
	}
	return (double) (uint16_t)(MaxScore - (uint16_t) 0x8000);
}

/*double swps3_alignmentShort2SSE(ProfileShort * query, const char * db,
		int dbLen, Options * options) {
	return 0;
}*/

EXPORT double swps3_alignmentShortSSE(ProfileShort * query, const char * db,
		int dbLen, Options * options) {

	/**********************************************************************
	 * This version of the code implements the idea presented in
	 *
	 ***********************************************************************
	 * Striped Smith-Waterman speeds database searches six times over other
	 * SIMD implementations
	 *
	 * Michael Farrar, Bioinformatics, 23(2), pp. 156-161, 2007
	 **********************************************************************/

	int i, j;
	uint16_t MaxScore = 0x8000;
	int segLength = (query->len + 7) / 8; /* the segment length */

	__m128i * loadOpt = query->loadOpt;
	__m128i * storeOpt = query->storeOpt;
	__m128i * rD = query->rD;
	__m128i * current_profile;
	__m128i * swap;

	__m128i vMinimums = _mm_set1_epi16(0x8000);

	__m128i vDelIncr = _mm_set1_epi16(options->gapExt);
	__m128i vDelFixed = _mm_set1_epi16(options->gapOpen);
	__m128i vBias = _mm_set1_epi16(query->bias);

	__m128i vMaxScore = vMinimums; /* vMaxScore = [0,0] */

	__m128i vProfile = vMinimums; /* the score profile */
	__m128i vStoreOpt; /* the new optimal score */
	__m128i vRD; /* the new row deletion score */
	__m128i vCD; /* the column deletion score */
	__m128i vTmp;

	if (options->gapExt <= options->gapOpen) {
		return swps3_alignmentShortSSE_lin(query, db, dbLen, options);
	}

	/* initialize the other arrays used for the dynProg code */
	/*********************************************************/
	for (i = 0; LIKELY(i < segLength); i++) {
		_mm_store_si128(loadOpt + i, vMinimums);
		_mm_store_si128(storeOpt + i, vMinimums);
		_mm_store_si128(rD + i, vMinimums);
	}

	/* looping through all the columns */
	/***********************************/
	for (j = 0; LIKELY(j < dbLen); j++) {

		/* compute the opt and cd score depending on the previous column */
		/*******************************************************************/
		/* set the column deletion score to zero, has to be fixed later on */
		vCD = vMinimums;

		/* set the opt score to the elements computed in the previous column */
		/* set the low of storeOpt to MaxS[j] */
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
		vStoreOpt = _mm_slli_si128(vStoreOpt, 2);
		vStoreOpt = _mm_insert_epi16(vStoreOpt, (int) 0x8000, 0);

		/* compute the current profile, depending on the character in s2 */
		/*****************************************************************/

		current_profile = query->profile + db[j] * segLength;
		/* swap the old optimal score with the new one */
		/***********************************************/
		swap = storeOpt;
		storeOpt = loadOpt;
		loadOpt = swap;

		/* main loop computing the max, precomputing etc. */
		/**************************************************/
		for (i = 0; LIKELY(i < segLength); i++) {

			vRD = _mm_load_si128(rD + i);
			vRD = _mm_adds_epi16(vRD, vDelIncr);
			vTmp = _mm_load_si128(loadOpt + i);
			vTmp = _mm_adds_epi16(vTmp, vDelFixed);
			vRD = _mm_max_epi16(vRD, vTmp);
			_mm_store_si128(rD + i, vRD);

			/* load the current profile */
			/*vProfile = _mm_movpi64_epi64(current_profile[i]);*/
			/*vProfile = _mm_loadl_epi64((__m128i*)(current_profile+i));*/
			/*#if (defined _WIN32 || defined __WIN32__)*/
			vProfile = _mm_load_si128(current_profile + i);
			/*#else
			 __asm__("MOVDQA (%1),%0" : "=x" (vProfile) : "r" (current_profile+i));
			 #endif*/
			/*vProfile = _mm_unpacklo_epi16(vProfile, _mm_xor_si128(vProfile,vProfile));*/
			vProfile = _mm_subs_epi16(vProfile, vBias);

			/* add the profile the prev. opt */
			vStoreOpt = _mm_adds_epi16(vStoreOpt, vProfile);

			/* update the maxscore found so far */
			vMaxScore = _mm_max_epi16(vMaxScore, vStoreOpt);

			/* compute the correct opt score of the cell */
			vStoreOpt = _mm_max_epi16(vStoreOpt, vCD);
			vStoreOpt = _mm_max_epi16(vStoreOpt, vRD);

			/* store the opt score of the cell */
			_mm_store_si128(storeOpt + i, vStoreOpt);

			/* precompute cd for next iteration */
			vStoreOpt = _mm_adds_epi16(vStoreOpt, vDelFixed);
			vCD = _mm_adds_epi16(vCD, vDelIncr);
			vCD = _mm_max_epi16(vCD, vStoreOpt);

			/* load precomputed opt for next iteration */
			vStoreOpt = _mm_load_si128(loadOpt + i);
		}

		for (i = 0; LIKELY(i < 8); ++i) {
			int k;
			/* compute the gap extend penalty for the current cell */
			vCD = _mm_slli_si128(vCD, 2);
			vCD = _mm_insert_epi16(vCD, 0x8000, 0);

			for (k = 0; LIKELY(k < segLength); ++k) {
				/* compute the current optimal value of the cell */
				vStoreOpt = _mm_load_si128(storeOpt + k);
				vStoreOpt = _mm_max_epi16(vStoreOpt, vCD);
				_mm_store_si128(storeOpt + k, vStoreOpt);

				/* precompute the scores for the next cell */
				vStoreOpt = _mm_adds_epi16(vStoreOpt, vDelFixed);
				vCD = _mm_adds_epi16(vCD, vDelIncr);

				if (UNLIKELY(
						!_mm_movemask_epi8(_mm_cmpgt_epi16(vCD, vStoreOpt))))
					goto shortcut;
			}
		}
		shortcut:
		/* store the new MaxScore for the next line block */
		/**************************************************/

		/* store the element of storeOpt in MaxS */
		vStoreOpt = _mm_load_si128(storeOpt + segLength - 1);
	}
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 8));
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 4));
	vMaxScore = _mm_max_epi16(vMaxScore, _mm_srli_si128(vMaxScore, 2));
	MaxScore = _mm_extract_epi16(vMaxScore, 0);
	if (MaxScore == 0x7fff) {
		return DBL_MAX;
	}
	return (double) (uint16_t)(MaxScore - (uint16_t) 0x8000);
}

EXPORT void swps3_freeProfileShortSSE(ProfileShort * profile) {
	free(profile);
}

