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

!define MULTIUSER_EXECUTIONLEVEL Highest
!define MULTIUSER_INSTALLMODE_DEFAULT_CURRENTUSER
!define MULTIUSER_MUI
!define MULTIUSER_INSTALLMODE_COMMANDLINE
!define MULTIUSER_INSTALLMODE_INSTDIR "SimNIBS-{{ version }}"
!define MULTIUSER_INSTALLMODE_FUNCTION correct_prog_files
!include MultiUser.nsh
!include FileFunc.nsh

; Modern UI installer stuff
!include "MUI2.nsh"
!define MUI_ABORTWARNING
!define MUI_ICON "gui_icon.ico"
!define MUI_UNICON "gui_icon.ico"

; UI pages
!insertmacro MUI_PAGE_WELCOME
!insertmacro MULTIUSER_PAGE_INSTALLMODE
!insertmacro MUI_PAGE_DIRECTORY
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


Name "SimNIBS ${PRODUCT_VERSION}"
OutFile "${INSTALLER_NAME}"
ShowInstDetails show

Var cmdLineInstallDir

Section -SETTINGS
  SetOutPath "$INSTDIR"
  SetOverwrite ifnewer
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

  ; Marker file for per-user install
  StrCmp $MultiUser.InstallMode CurrentUser 0 +3
    FileOpen $0 "$INSTDIR\${USER_INSTALL_MARKER}" w
    FileClose $0
    SetFileAttributes "$INSTDIR\${USER_INSTALL_MARKER}" HIDDEN

  ; Run Scripts
  ; These steps rely on the fix_entrypoints.py and postinstall_simnibs.py being moved to the simnibs_env dir
  ; The sitecustomize.py also needs to be moved to simnibs_env/Lib/site-packages in order for the postinstall_simnibs.py to run
  DetailPrint "Fixing Scripts"
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
  nsExec::ExecToLog '"$INSTDIR\simnibs_env\Scripts\postinstall_simnibs.exe" -d "$INSTDIR" --copy-matlab --setup-links --no-extra-coils'
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
  WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayName" "${PRODUCT_NAME}"
  WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "InstallLocation" "$INSTDIR"
  WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayIcon" "$INSTDIR\${PRODUCT_ICON}"
  WriteRegStr SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "DisplayVersion" "${PRODUCT_VERSION}"
  WriteRegDWORD SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "NoModify" 1
  WriteRegDWORD SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}" \
                   "NoRepair" 1

  ; Check if we need to reboot
  IfRebootFlag 0 noreboot
    MessageBox MB_YESNO "A reboot is required to finish the installation. Do you wish to reboot now?" \
                /SD IDNO IDNO noreboot
      Reboot
  noreboot:

  StrCpy $0 "$INSTDIR\install.log"
  Push $0
  Call DumpLog

SectionEnd

Section "Uninstall"
  SetRegView 64
  SetShellVarContext all
  IfFileExists "$INSTDIR\${USER_INSTALL_MARKER}" 0 +3
    SetShellVarContext current
    Delete "$INSTDIR\${USER_INSTALL_MARKER}"

  ; Run the postinstall uninstaller to remove from PATH and Start Menu
  ; The uninstall_simnibs.cmd is created by the postinstall script
  nsExec::ExecToLog '"$INSTDIR\uninstall_simnibs.cmd"'

  Delete "$INSTDIR\${PRODUCT_ICON}"
  ; Uninstall directories
  RMDir /r "$INSTDIR\simnibs_env"
  RMDir /r "$INSTDIR\documentation"
  Delete "$INSTDIR\install.log"
  Delete "$INSTDIR\uninstall.exe"
  Delete "$INSTDIR\uninstall_simnibs.cmd"

  RMDir $INSTDIR
  DeleteRegKey SHCTX "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
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
  StrCpy $INSTDIR "$PROFILE\${PRODUCT_NAME}"
FunctionEnd

Function un.onInit
  !insertmacro MULTIUSER_UNINIT
FunctionEnd

Function correct_prog_files
  ; The multiuser machinery doesn't know about the different Program files
  ; folder for 64-bit applications. Override the install dir it set.
  StrCmp $MultiUser.InstallMode AllUsers 0 +2
    StrCpy $INSTDIR "$PROGRAMFILES64\${MULTIUSER_INSTALLMODE_INSTDIR}"
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