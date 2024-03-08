.. _install_throubleshooting:

Troubleshooting
===============

* Antivirus software and other security features might interfere with the installation.

* if you are seeing the error messages such as

    * :code:`.plugin: Could not find the Qt platform plugin "xcb" in ""`
    * :code:`.plugin: Could not find the Qt platform plugin "windows" in ""`

    You likely have `non-ASCII charactes <https://en.wikipedia.org/wiki/ASCII>`_ (e.g. Chinese characters) in your installation path. Please change the installation path so that it only includes ASCII characters.

* Often times, :ref:`installing with the conda package manager <conda-install>` can help overcome problems with the installer. Please try it out

* **glError** messages when opening the GUI can often be solved by updating video drivers.

* If you are stuck, please :ref:`contact us <contact>` and attach the :file:`simnibs_install_log.txt` file, which can be found in your SimNIBS installation directory
