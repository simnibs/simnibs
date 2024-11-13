## Copyright (C) 2015 Markus Bergholz <markuman@gmail.com>
## Co-author  Nicholas R. Jankowski <jankowskin@asme.org> **
##   ** U.S. government employee, contributions are public domain
##      with no personal assertion of copyright
## prepareIdx routine Copyright (C) 2015 Peter John Acklam <pjacklam@gmail.com>, 
## used with permission.
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {} repelem (@var{X}, @var{R})
## @deftypefnx {} {} repelem (@var{X}, @var{R1}, @dots{}, @var{Rn})
## Construct an array of repeated elements from X.
##
## @var{X} must be a scalar, a vector or a n-dimensional array. 
##
## @var{Rn} must either be scalar or a vectror with the same number of
## elements as the size of dimension n of the array to be replicated. @var{Rn} may
## be either row or column vectors.
##
## repelem may only be called with a single @var{R} if @var{X} is a scalar or 
## vector. @var{X} may be either a row or column vector, and repelem will return
## a replicated vector with the same orientation.  Non-vector @var{X}'s require
## at least two @var{Rn}'s be specified. (Using repelem with a vector @var{X} and 
## a vector for @var{Rn} is equivalent to Run Length Decoding.)
##
##@example
##A = [1 2 3 4 5];
##B = [2 1 0 1 2];
##repelem(A,B)
##  ==>    1   1   2   4   5   5
##@end example
##@example
##A = magic(3)
##  ==>    8   1   6
##         3   5   7
##         4   9   2
##B1 = [1 2 3];
##B2 = 2;
##repelem(A,B1,B2)
##  ==>    8   8   1   1   6   6
##         3   3   5   5   7   7
##         3   3   5   5   7   7
##         4   4   9   9   2   2
##         4   4   9   9   2   2
##         4   4   9   9   2   2
##@end example
##
## More @var{Rn} may be specified than the number of dimensions of @var{X}.
## Those excess @var{Rn} must be scalars (because @var{X}'s size in those 
## dimensions is only 1), and @var{X} will be replicated in those dimension 
## accordingly. 
##  
##@example
## A = [1 2 3 4 5];
## B1 = 2;
## B2 = [2 1 3 0 2];
## B3 = 3;
## repelem(A,B1,B2,B3)
##  ==>    ans(:,:,1) =
##            1   1   2   3   3   3   5   5
##            1   1   2   3   3   3   5   5
##         
##         ans(:,:,2) =
##         
##            1   1   2   3   3   3   5   5
##            1   1   2   3   3   3   5   5
##         
##         ans(:,:,3) =
##            1   1   2   3   3   3   5   5
##            1   1   2   3   3   3   5   5
##@end example
##
## @var{Rn} must be specified in order. A placeholder of 1 can be used in 
## deminsions not needing replication. [] inputs will produce an error. 
##
##@example
## repelem([-1 0;0 1],1,2,1,2)
##  ==>  ans(:,:,1,1) =
##         -1  -1   0   0
##          0   0   1   1
##
##       ans(:,:,1,2) =
##         -1  -1   0   0
##          0   0   1   1 
##@end example
##
## If fewer @var{Rn} are given than dimensions in @var{X}, repelem will assume a 
## @var{Rn} of 1 for those dimensions.
##
##@example
## A = cat(3,[-1 0;0 1],[-1 0;0 1])
##  ==>  ans(:,:,1) =
##         -1   0
##          0   1
##
##       ans(:,:,2) =
##         -1   0
##          0   1
##
## repelem(A,2,3)
##  ==>  ans(:,:,1) =
##         -1  -1  -1   0   0   0
##         -1  -1  -1   0   0   0
##          0   0   0   1   1   1
##          0   0   0   1   1   1
##
##       ans(:,:,2) =
##         -1  -1  -1   0   0   0
##         -1  -1  -1   0   0   0
##          0   0   0   1   1   1
##          0   0   0   1   1   1
##@end example 
##
## repelem will preserve the class of @var{X}, and works with strings, NA and
## NAN inputs. Any @var{Rn} = 0 will produce an empty array. 
##
##@example
## repelem('Octave',2,3)
##  ==>    OOOccctttaaavvveee
##         OOOccctttaaavvveee
##
##repelem([1 2 3;1 2 3],2,0)
##  ==>    [](4x0)
##@end example
##
## @seealso{cat korn repmat}
## @end deftypefn


function retval = repelem(element, varargin)

  if (nargin <= 1)
    error("repelem: Not enough input arguments")
    print_usage ();
    
  elseif (nargin == 2)
  
    v = varargin{1};

    if (isscalar (v))

      if (iscolumn (element))
        # element values repeated v times in a col vector
        retval = element (:, ones (v, 1))'(:);

      elseif (isrow (element))
        # element values repeated v times in a row vector
        retval = element (ones (v, 1), :)(:)';

      else
        error (["repelem: %gD Array objects require %g or more input " ...
               "arguments, only %g given"], ...
                ndims (element), ndims (element) + 1, nargin);
        print_usage ();

      endif

    elseif (isvector (element)) && (isvector (v))

      if (numel (v) == numel (element))
        # vector element with vector varargin. basic run-length decoding in
        # function prepareIdx returned idx2 as a row vect. of element indices
        # in right position
        idx2 = prepareIdx (v);
        # fills with element values, direction matches element. 
        retval  = element (idx2);

      else # catch unequal element counts
        error (["repelem: varargin must either be scalar or have the same" ...
                " number of elements as the vector to be replicated"]);
        print_usage ();

      endif
      
    else # catch any arrays passed to element or varargin with nargin==2
      error (["repelem: when called with only two inputs they must be" ...
              " either scalars or vectors."]);
      print_usage ();

    endif
  
  
  elseif (nargin == 3)  #can simplify for known dimension count
    
    # avoid repeated function calls
    elsize = size (element);
    # 'numel' or 'length' faster than isvector in cellfun
    scalarv = (cellfun ('numel', varargin) == 1);
    nonscalarv = ~scalarv;
    
    ##INPUT CHECK

    #1:check that all varargin are either scalars or vectors, no arrays.
    # isvector gives true for scalars.
    # (Faster here with only two to avoid cellfun)
    if (~ (isvector (varargin{1}) && (isvector (varargin{2}))))
      error ("repelem: varargin must be all scalars or vectors");
      print_usage ();

    #2: check that the ones that are vectors have the right length.
    elseif (any (~ (cellfun ('numel', varargin (nonscalarv)) == elsize (nonscalarv))))
      error (["repelem: varargin(n) must either be scalar or have the same" ...
              " number of elements as the size of dimension n of the array" ...
              " to be replicated"]);
      print_usage ();

    endif 
    
    #Create index arrays to pass to element 
    ##(no slower passing to prepareIdx than checking and doing scalars directly)
    idx1 = prepareIdx (varargin{1}, elsize (1));
    idx2 = prepareIdx (varargin{2}, elsize (2));
  
    #the : at the end takes care of size(element)>2
    retval = element (idx1, idx2, :);

  else  #if (nargin > 3) **no need for elseif
  
    # avoid repeated function calls
    elsize = size (element);
    eldims = numel (elsize);
    vasize = nargin - 1; #numel(varargin);
    %maxDim = max(eldims,vasize);
    dimsWithBoth = min (eldims, vasize);
    nonscalarv   = ~cellfun (@isscalar, varargin);

    ## INPUT CHECK

    # 1: that they are all scalars or vectors. isvector gives true for scalars.
    if (~all (cellfun (@isvector, varargin)))
      error ("repelem: varargin must be all be scalars or vectors");
      print_usage ();

    # 2: catch any vectors thrown at trailing singletons, which should only have scalars
    elseif (max (find (nonscalarv)) > eldims)
      error ("repelem: varargin(n) for trailing singleton dimensions must be scalar");
      print_usage ();

    # 3: that the ones that are vectors have the right length.
    elseif (any (~ (cellfun ('numel', varargin (nonscalarv)) == elsize (nonscalarv))))
      error (["repelem: varargin(n) must either be scalar or have the same" ...
              " number of elements as the size of dimension n of the array" ...
              " to be replicated"]);
      print_usage ();

    endif 
    
    # first, preallocate idx which will contain index array to be put into element
    idx = cell (1, vasize);

    # use prepareIdx() to fill indices for each dimension that could be a scalar or vector
    idx (1:dimsWithBoth) = cellfun (@prepareIdx, varargin (1:dimsWithBoth), ...
                           num2cell (elsize (1:dimsWithBoth)), ...
                           'UniformOutput', false);

    # if there are more varargin inputs than element dimensions, input tests have 
    # verified they are just scalars, so add [1 1 1 1 1... 1] to those dims to 
    # perform concatenation in those dims.
    
    ###can cellfun check speed against for loop.
    if (vasize > eldims)
      idx (eldims + 1:vasize) = cellfun ('ones', {1}, ...
                                {varargin(eldims + 1:end){:}}, ...
                                'UniformOutput', false);
    endif
    
    # use completed idx to specify repitition of element values in all dimensions
    # trailing : will take care of any case where eldims > vasize.
    retval = element (idx{:}, :);
  
  endif

endfunction


function idx = prepareIdx (v, elsize_n)
# returns a row vector of indices prepared for replicating.

  if (isscalar (v))
    # will always return row vector
    idx = [1:elsize_n](ones (v, 1), :)(:)';

  else
  
    # works for row or column vector. idx2 output will be a row vector.
    # gets ending position for each element item
    idx_temp = cumsum (v);
    # row vector with enough space for output
    idx (1:idx_temp (end)) = 0;
    # sets starting position of each element to 1
    idx (idx_temp (1:end - 1) + 1) = 1;
    # sets starting position of each element to 1
    idx(1) = 1;
    # with prepared index
    idx = (find (v ~= 0))(cumsum (idx));

  endif
  
endfunction

%%%tests for help examples
%!assert (repelem([1 2 3 4 5], [2 1 0 1 2]), [1 1 2 4 5 5])
%!assert (repelem(magic(3), [1 2 3],2), ...
%!  [8 8 1 1 6 6;3 3 5 5 7 7;3 3 5 5 7 7;4 4 9 9 2 2;4 4 9 9 2 2;4 4 9 9 2 2])
%!assert (repelem([1 2 3 4 5],2,[2 1 3 0 2],3),repmat([1 1 2 3 3 3 5 5],2,1,3))
%!assert (repelem([-1 0;0 1],1,2,1,2), repmat([-1 -1 0 0; 0 0 1 1],1,1,1,2))
%!assert (repelem(cat(3,[-1 0 ; 0 1],[-1 0 ; 0 1]),2,3), ...
%!  repmat([-1 -1 -1 0 0 0;-1 -1 -1 0 0 0;0 0 0 1 1 1;0 0 0 1 1 1],1,1,2))
%!assert (repelem('Octave', 2,3), ['OOOccctttaaavvveee';'OOOccctttaaavvveee'])
%
%%%nargin ==2 tests
%!assert (repelem([-1 0 1], 2), [-1 -1 0 0 1 1])
%!assert (repelem([-1 0 1]', 2), [-1; -1; 0; 0; 1; 1;])
%!assert (repelem([-1 0 1], [1 2 1]), [-1 0 0 1]) 
%!assert (repelem([-1 0 1]', [1 2 1]), [-1; 0; 0; 1])
%!assert (repelem([1 2 3 4 5]', [2 1 0 1 2]), [1 1 2 4 5 5]')
%
%%%nargin ==3 tests
%!assert (repelem([1 0;0 -1], 2, 3), ...
%!  [1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 -1 -1 -1;0 0 0 -1 -1 -1])
%!assert (repelem([1 0; 0 -1], 1,[3 2]), [1 1 1 0 0;0 0 0 -1 -1])
%!assert (repelem([1 0; 0 -1], 2,[3 2]), [1 1 1 0 0;1 1 1 0 0;0 0 0 -1 -1;0 0 0 -1 -1])
%!assert (repelem(cat(3,[1 0; 0 -1],[1 0;0 -1]), 1,[3 2]), ...
%!  repmat([1 1 1 0 0 ; 0 0 0 -1 -1],1,1,2))
%!assert (repelem([1 0; 0 -1], [3 2], 1), [1 0;1 0;1 0;0 -1;0 -1])
%!assert (repelem([1 0; 0 -1], [3 2], 2), ...
%!  [1 1 0 0;1 1 0 0;1 1 0 0;0 0 -1 -1;0 0 -1 -1])
%!assert (repelem([1 0; 0 -1], [2 3] ,[3 2]), ...
%!  [1 1 1 0 0;1 1 1 0 0;0 0 0 -1 -1;0 0 0 -1 -1;0 0 0 -1 -1])
%!assert (repelem(cat(3,[1 1 1 0;0 1 0 0],[1 1 1 1;0 0 0 1],[1 0 0 1;1 1 0 1]), 2, 3), ...
%!  cat(3,[1 1 1 1 1 1 1 1 1 0 0 0;1 1 1 1 1 1 1 1 1 0 0 0;0 0 0 1 1 1 0 0 0 0 0 0;0 0 0 1 1 1 0 0 0 0 0 0], ...
%!  [1 1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0 0 1 1 1;0 0 0 0 0 0 0 0 0 1 1 1], ...
%!  [1 1 1 0 0 0 0 0 0 1 1 1;1 1 1 0 0 0 0 0 0 1 1 1;1 1 1 1 1 1 0 0 0 1 1 1;1 1 1 1 1 1 0 0 0 1 1 1]))
%!assert (repelem(cat(3,[1 1 1 0;0 1 0 0],[1 1 1 1;0 0 0 1],[1 0 0 1;1 1 0 1]), ...
%!  2, [3 3 3 3]), ...
%!  cat(3,[1 1 1 1 1 1 1 1 1 0 0 0;1 1 1 1 1 1 1 1 1 0 0 0;0 0 0 1 1 1 0 0 0 0 0 0;0 0 0 1 1 1 0 0 0 0 0 0], ...
%!  [1 1 1 1 1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0 0 1 1 1;0 0 0 0 0 0 0 0 0 1 1 1], ...
%!  [1 1 1 0 0 0 0 0 0 1 1 1;1 1 1 0 0 0 0 0 0 1 1 1;1 1 1 1 1 1 0 0 0 1 1 1;1 1 1 1 1 1 0 0 0 1 1 1]))
%!assert (repelem([1 2 3 4 5], 2,[2 1 2 0 2]), [1 1 2 3 3 5 5;1 1 2 3 3 5 5])
%
%%%nargin > 3 tests
%!assert (repelem([1 0;0 -1], 2, 3, 4), ...
%!  cat(3,[1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 -1 -1 -1;0 0 0 -1 -1 -1], ...
%!  [1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 -1 -1 -1;0 0 0 -1 -1 -1], ...
%!  [1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 -1 -1 -1;0 0 0 -1 -1 -1], ...
%!  [1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 -1 -1 -1;0 0 0 -1 -1 -1]))
%!assert (repelem(repmat([-1 0;0 1],1,1,2,3),2,2,2), ...
%!  repmat([-1 -1 0 0;-1 -1 0 0;0 0 1 1; 0 0 1 1],1,1,4,3))
%!assert (repelem(repmat([-1 0;0 1],1,1,2,3),[2 2],[2 2],2), ...
%!  repmat([-1 -1 0 0;-1 -1 0 0;0 0 1 1; 0 0 1 1],1,1,4,3))
%!assert (repelem(repmat([-1 0;0 1],1,1,2,3),2,2,2,2,2), ...
%!  repmat([-1 -1 0 0;-1 -1 0 0;0 0 1 1; 0 0 1 1],1,1,4,6,2))
%!assert (repelem([1,0,-1;-1,0,1],[2 3],[2 3 4],2), ...
%!  cat(3,[1 1 0 0 0 -1 -1 -1 -1;1 1 0 0 0 -1 -1 -1 -1;-1 -1 0 0 0 1 1 1 1;-1 -1 0 0 0 1 1 1 1;-1 -1 0 0 0 1 1 1 1], ...
%!  [1 1 0 0 0 -1 -1 -1 -1;1 1 0 0 0 -1 -1 -1 -1;-1 -1 0 0 0 1 1 1 1;-1 -1 0 0 0 1 1 1 1;-1 -1 0 0 0 1 1 1 1]))
%!assert (repelem([1 2 3;4 5 6],[0 2],2,2), repmat([4 4 5 5 6 6],2,1,2))
%
%%%nargin <= 1 error tests
%!error  (repelem())
%!error  (repelem([]))
%!error  (repelem(5))
%!error  (repelem(5,[]))
%!error  (repelem([1 2 3 3 2 1]))
%!error  (repelem([1 2 3; 3 2 1]))
%
%%%nargin == 2 error tests
%!error  (repelem([1 2 3; 3 2 1],[]))
%!error  (repelem([1 2 3; 3 2 1],2))
%!error  (repelem([1 2 3; 3 2 1],2))
%!error  (repelem([1 2 3; 3 2 1],[1 2 3]))
%!error  (repelem([1 2 3; 3 2 1],[1 2 3]'))
%!error  (repelem([1 2 3; 3 2 1],[1 2 2 1]))
%!error  (repelem([1 2 3; 3 2 1],[1 2 3;4 5 6]))
%!error  (repelem([1 2 3 4 5],[1 2 3 4 5;1 2 3 4 5]))
%
%%%nargin == 3 error tests
%!error  (repelem([1 2 3; 3 2 1], 1, [1 2;1 2]))
%!error  (repelem([1 2 3; 3 2 1], 1, [1 2]))
%!error  (repelem([1 2 3; 3 2 1], 2, []))
%!error  (repelem([1 2 3; 3 2 1], [1 2 3], [1 2 3]))
%!error  (repelem([1 2 3; 3 2 1], [1 2 3], [1 2 3 4]))
%!error  (repelem([1 2 3; 3 2 1], [1 2], [1 2 3 4]))
%
%%%nargin > 3 error tests
%!error  (repelem([1 2 3; 3 2 1], 1, [1 2;1 2],1,2,3))
%!error  (repelem([1 2 3; 3 2 1], [],1,2,3))
%!error  (repelem([1 2 3; 3 2 1], [1 2], [1 2 3],1,2,[1 2;1 2]))
%!error  (repelem([1 2 3; 3 2 1], [1 2 3], [1 2 3],1,2))
%!error  (repelem([1 2 3; 3 2 1], [1 2], [1 2 3 4],1,2))