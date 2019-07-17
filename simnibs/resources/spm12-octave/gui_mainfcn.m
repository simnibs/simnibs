function varargout = gui_mainfcn (gui_State, varargin)

  if (numel (varargin) == 0)

    ## Open the figure
    copies = ifelse (gui_State.gui_Singleton, "reuse", "new");
    if (! isempty (gui_State.gui_LayoutFcn))
      F = feval (gui_State.gui_LayoutFcn, copies);
    else
      figname = which (gui_State.gui_Name);
      [p, n]  = fileparts (figname);
      figname = fullfile (p, n);
      F       = openfig (figname, copies);
    endif

    ## Store handles in guidata
    D = guidata (F);
    H = guihandles (F);
    guidata (F, H);

    ## Call gui_OpeningFcn
    if (! isempty (gui_State.gui_OpeningFcn))
      feval (gui_State.gui_OpeningFcn, F, [], H);
    endif

    ## Call gui_OutputFcn
    H = guidata (F);
    if (nargout > 0)
      [varargout{1:nargout}] = feval (gui_State.gui_OutputFcn, F, [], H);
    else
      feval (gui_State.gui_OutputFcn, F, [], H);
    endif

  else

    ## Callbacks
    if (nargout > 0)
      [varargout{1:nargout}] = feval (gui_State.gui_Callback, varargin{2:end});
    else
      feval (gui_State.gui_Callback, varargin{2:end});
    endif

  endif

endfunction
