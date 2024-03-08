/* Generic x86 gmp-mparam.h -- Compiler/machine parameter header file.

Copyright 1991, 1993, 1994, 2000, 2001, 2002 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

/*
 * This gmp-mparam.h is a wrapper include file for the original gmp-mparam.h, 
 * which has been renamed to gmp-mparam-<arch>.h. There are conflicts for the
 * original gmp-mparam.h on multilib systems, which result from arch-specific
 * configuration options. Please do not use the arch-specific file directly.
 *
 * Copyright (C) 2006 Red Hat, Inc.
 * Thomas Woerner <twoerner@redhat.com>
 */

#ifdef gmp_mparam_wrapper_h
#error "gmp_mparam_wrapper_h should not be defined!"
#endif
#define gmp_mparam_wrapper_h

#if defined(__arm__)
#include "gmp-mparam-arm.h"
#elif defined(__i386__)
#include "gmp-mparam-i386.h"
#elif defined(__ia64__)
#include "gmp-mparam-ia64.h"
#elif defined(__powerpc64__)
#include "gmp-mparam-ppc64.h"
#elif defined(__powerpc__)
#include "gmp-mparam-ppc.h"
#elif defined(__s390x__)
#include "gmp-mparam-s390x.h"
#elif defined(__s390__)
#include "gmp-mparam-s390.h"
#elif defined(__x86_64__)
#include "gmp-mparam-x86_64.h"
#elif defined(__alpha__)
#include "gmp-mparam-alpha.h"
#elif defined(__sh__)
#include "gmp-mparam-sh.h"
#elif defined(__sparc__) && defined (__arch64__)
#include "gmp-mparam-sparc64.h"
#elif defined(__sparc__)                      
#include "gmp-mparam-sparc.h"
#else
#error "The gmp-devel package is not usable with the architecture."
#endif

#undef gmp_mparam_wrapper_h
