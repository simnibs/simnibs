!define PRODUCT_NAME "SimNIBS-{{ version }}"
!define PRODUCT_VERSION "{{ full_version }}"
!define INSTALLER_NAME "simnibs_installer_windows.exe"
!define PRODUCT_ICON "gui_icon.ico"

; Marker file to tell the uninstaller that it's a user installation
!define USER_INSTALL_MARKER _user_install_marker

SetCompressor lzma

!if "${NSIS_PACKEDVERSION}" >= 0x03000000
  Unicode true
  ManifestDPIAware true
!endif

RequestExecutionLevel user
!include FileFunc.nsh

; Modern UI installer stuff
!include "MUI2.nsh"
!define MUI_ABORTWARNING
!define MUI_ICON "gui_icon.ico"
!define MUI_UNICON "gui_icon.ico"

; UI pages
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
Page Custom ShowOptions
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_LANGUAGE "English"

; Logging
!define /IfNDef LVM_GETITEMCOUNT 0x1004
!define /IfNDef LVM_GETITEMTEXTA 0x102D
!define /IfNDef LVM_GETITEMTEXTW 0x1073
!if "${NSIS_CHAR_SIZE}" > 1
!define /IfNDef LVM_GETITEMTEXT ${LVM_GETITEMTEXTW}
!else
!define /IfNDef LVM_GETITEMTEXT ${LVM_GETITEMTEXTA}
!endif

!include LogicLib.nsh
!include x64.nsh

var hwnd_dialog
var hwnd_gmshoptions
var hwnd_systempath
var hwnd_shortcut
var hwnd_associatefiles
var hwnd_messagestring
var messagestring
var postinstall_options
var python_note

Function ShowOptions
    !insertmacro MUI_HEADER_TEXT "Post-Installation Options" ""
	nsDialogs::Create 1018
	Pop $hwnd_dialog
	
	${NSD_CreateCheckbox} 0 0 50% 6% "Copy gmsh_options file"
		Pop $hwnd_gmshoptions
		${NSD_Check} $hwnd_gmshoptions
		${NSD_OnClick} $hwnd_gmshoptions UpdateText
		
	${NSD_CreateCheckbox} 0 10% 50% 6% "Add SimNIBS to the system PATH"
        Pop $hwnd_systempath
        ${NSD_Check} $hwnd_systempath
        ${NSD_OnClick} $hwnd_systempath UpdateText
        
   	${NSD_CreateCheckbox} 0 20% 50% 6% "Add SimNIBS shortcut icons"
        Pop $hwnd_shortcut
        ${NSD_Check} $hwnd_shortcut
        ${NSD_OnClick} $hwnd_shortcut UpdateText
 
    ${NSD_CreateCheckbox} 0 30% 50% 6% "Associate .msh, .stl and .geo files"
        Pop $hwnd_associatefiles
        ${NSD_Check} $hwnd_associatefiles
        ${NSD_OnClick} $hwnd_associatefiles UpdateText
        
    ${NSD_CreateLabel} 0 50% 100% 30% ""
        Pop $hwnd_messagestring
     
    StrCpy $python_note "NOTE: The installer will run SimNIBS configuration scripts. Windows might ask therefore ask"
    StrCpy $python_note "$python_note for your OK to run python during installation."
    ${NSD_CreateLabel} 0 80% 100% 20% $python_note
    
	nsDialogs::Show
FunctionEnd

Function UpdateText
    StrCpy $messagestring ""
    StrCpy $postinstall_options ""
    
    ${NSD_GetState} $hwnd_gmshoptions $0
    ${If} $0 == 0
        StrCpy $messagestring "$messagestring Visualization of head models and simulation results will be affected!$\n"
        StrCpy $postinstall_options "$postinstall_options --no-copy-gmsh-options"
    ${EndIf}
    
    ${NSD_GetState} $hwnd_systempath $0
    ${If} $0 == 0
        StrCpy $messagestring "$messagestring SimNIBS functions will not be callable from the terminal!$\n"
        StrCpy $postinstall_options "$postinstall_options --no-add-to-path"
    ${EndIf}
    
    ${NSD_GetState} $hwnd_shortcut $0
    ${If} $0 == 0
        StrCpy $messagestring "$messagestring SimNIBS icons will not be found in the Start Menu!$\n"
        StrCpy $postinstall_options "$postinstall_options --no-shortcut-icons"
    ${EndIf}
    
    ${NSD_GetState} $hwnd_associatefiles $0
    ${If} $0 == 0
        StrCpy $messagestring "$messagestring Models and simulation results will not be associated with Gmsh!$\n"
        StrCpy $postinstall_options "$postinstall_options --no-associate-files"
    ${EndIf}
    
    ${NSD_SetText} $hwnd_messagestring $messagestring
FunctionEnd


Name "SimNIBS ${PRODUCT_VERSION}"
OutFile "${INSTALLER_NAME}"
ShowInstDetails show

Section -SETTINGS
  SetOutPath "$INSTDIR"
  SetOverwrite on
SectionEnd


Section "!${PRODUCT_NAME}" sec_app
  SetRegView 64
  ;Check byteness
  ${IfNot} ${RunningX64}
    MessageBox MB_ICONSTOP "32 Bit Windows Detected. Can not install SimNIBS"
    Abort
  ${EndIf} 
  SectionIn RO
  File ${PRODUCT_ICON}
  SetOutPath "$INSTDIR\simnibs_env"
    File /r "simnibs_env\*.*"
  SetOutPath "$INSTDIR"

  SetOutPath "$INSTDIR\documentation"
    File /r "documentation\*.*"
  SetOutPath "$INSTDIR"

  ; Run Scripts
  ; These steps rely on the fix_entrypoints.py and postinstall_simnibs.py being moved to the simnibs_env dir
  ; The sitecustomize.py also needs to be moved to simnibs_env/Lib/site-packages in order for the postinstall_simnibs.py to run
  DetailPrint "Fixing Scripts"
  DetailPrint 'Running "$INSTDIR\simnibs_env\python.exe" "$INSTDIR\simnibs_env\fix_entrypoints.py" "$INSTDIR\simnibs_env\Scripts"  "$INSTDIR\simnibs_env"'
  nsExec::ExecToLog '"$INSTDIR\simnibs_env\python.exe" "$INSTDIR\simnibs_env\fix_entrypoints.py" "$INSTDIR\simnibs_env\Scripts"  "$INSTDIR\simnibs_env"'
  Pop $0
  ${IfNot} $0 == 0
      MessageBox MB_ICONSTOP "There was an error installing SimNIBS"
      StrCpy $0 "$INSTDIR\install.log"
      Push $0
      Call DumpLog
      Abort
  ${EndIf}
  
  DetailPrint "Installing the SimNIBS package..."
  ${If} ${Silent}
    StrCpy $postinstall_options ""
  ${EndIf}
  DetailPrint 'Running "$INSTDIR\simnibs_env\Scripts\postinstall_simnibs.exe" --silent --force -d "$INSTDIR" --setup-links $postinstall_options'
  nsExec::ExecToLog '"$INSTDIR\simnibs_env\Scripts\postinstall_simnibs.exe" --silent --force -d "$INSTDIR" --setup-links $postinstall_options'
  Pop $0
  ${IfNot} $0 == 0
      MessageBox MB_ICONSTOP "There was an error installing SimNIBS"
      StrCpy $0 "$INSTDIR\install.log"
      Push $0
      Call DumpLog
      Abort
  ${EndIf}


  WriteUninstaller $INSTDIR\uninstall.exe
  ; Add ourselves to Add/remove programs
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayName" "${PRODUCT_NAME}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "InstallLocation" "$INSTDIR"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayIcon" "$INSTDIR\${PRODUCT_ICON}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayVersion" "${PRODUCT_VERSION}"
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "NoRepair" 1

  StrCpy $0 "$INSTDIR\install.log"
  Push $0
  Call DumpLog

SectionEnd

Section "Uninstall"
  SetRegView 64
  ; Run the postinstall uninstaller to remove from PATH and Start Menu
  ; The uninstall_simnibs.cmd is created by the postinstall script
  nsExec::ExecToLog '"$INSTDIR\uninstall_simnibs.cmd" --silent'

  Delete "$INSTDIR\${PRODUCT_ICON}"
  ; Uninstall directories
  RMDir /r "$INSTDIR\simnibs_env"
  RMDir /r "$INSTDIR\documentation"
  Delete "$INSTDIR\install.log"
  Delete "$INSTDIR\uninstall.exe"
  Delete "$INSTDIR\uninstall_simnibs.cmd"
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"

  RMDir $INSTDIR
SectionEnd


; Functions

Function .onMouseOverSection
    ; Find which section the mouse is over, and set the corresponding description.
    FindWindow $R0 "#32770" "" $HWNDPARENT
    GetDlgItem $R0 $R0 1043 ; description item (must be added to the UI)

    StrCmp $0 ${sec_app} "" +2
      SendMessage $R0 ${WM_SETTEXT} 0 "STR:${PRODUCT_NAME}"

FunctionEnd

Function .onInit
  ; Change default to HOME folder
  InitPluginsDir
  InitPluginsDir
  ${If} $INSTDIR == ""
    StrCpy $INSTDIR "$PROFILE\${PRODUCT_NAME}"
  ${EndIf}
FunctionEnd

; Logging function from https://nsis.sourceforge.io/Dump_log_to_file
Function DumpLog
  Exch $5
  Push $0
  Push $1
  Push $2
  Push $3
  Push $4
  Push $6
  FindWindow $0 "#32770" "" $HWNDPARENT
  GetDlgItem $0 $0 1016
  StrCmp $0 0 exit
  FileOpen $5 $5 "w"
  StrCmp $5 "" exit
    SendMessage $0 ${LVM_GETITEMCOUNT} 0 0 $6
    System::Call '*(&t${NSIS_MAX_STRLEN})p.r3'
    StrCpy $2 0
    System::Call "*(i, i, i, i, i, p, i, i, i) i  (0, 0, 0, 0, 0, r3, ${NSIS_MAX_STRLEN}) .r1"
    loop: StrCmp $2 $6 done
      System::Call "User32::SendMessage(i, i, i, i) i ($0, ${LVM_GETITEMTEXT}, $2, r1)"
      System::Call "*$3(&t${NSIS_MAX_STRLEN} .r4)"
      !ifdef DumpLog_As_UTF16LE
      FileWriteUTF16LE ${DumpLog_As_UTF16LE} $5 "$4$\r$\n"
      !else
      FileWrite $5 "$4$\r$\n" ; Unicode will be translated to ANSI!
      !endif
      IntOp $2 $2 + 1
      Goto loop
    done:
      FileClose $5
      System::Free $1
      System::Free $3
  exit:
    Pop $6
    Pop $4
    Pop $3
    Pop $2
    Pop $1
    Pop $0
    Pop $5
FunctionEnd
