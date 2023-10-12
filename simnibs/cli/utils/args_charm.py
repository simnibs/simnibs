from .helpers import CommandLineArgument


subid = CommandLineArgument(
    ["subID"],
    dict(
        nargs="?",
        help="""Subject ID. The m2m_{sub_ID}
            folder will be created in the current working directory to
            store the result files.""",
    ),
)

primary_image = CommandLineArgument(["T1"], dict(nargs="?", help="T1-weighted image"))
secondary_image = CommandLineArgument(
    ["T2"], dict(nargs="?", help="T2-weighted image (optional)")
)
register_t2 = CommandLineArgument(
    ["--registerT2"],
    dict(
        action="store_true",
        default=False,
        help="Register T2- to T1-weighted image",
    ),
)
init_atlas = CommandLineArgument(
    ["--initatlas"],
    dict(
        action="store_true",
        default=False,
        help="""Affine registration of atlas to input images
                        (Note:T1-weighted image has to be supplied if no
                        T2-weighted image is used and --registerT2 is thus
                        skipped)""",
    ),
)
segment = CommandLineArgument(
    ["--segment"],
    dict(
        action="store_true",
        default=False,
        help="""Run segmentation to create label image,
                    reconstruct the middle cortical surfaces, and create
                    the registrations to the fsaverage and MNI templates""",
    ),
)
mesh = CommandLineArgument(
    ["--mesh"],
    dict(
        action="store_true",
        default=False,
        help="Create the head mesh from the label image",
    ),
)
surfaces = CommandLineArgument(
    ["--surfaces"],
    dict(
        action="store_true",
        default=False,
        help="Create central cortical surfaces from the label image",
    ),
)
forcerun = CommandLineArgument(
    ["--forcerun"],
    dict(
        action="store_true",
        default=False,
        help="""Overwrite existing m2m_{subID} folder instead
                    of throwing an error""",
    ),
)
skip_register_t2 = CommandLineArgument(
    ["--skipregisterT2"],
    dict(
        action="store_true",
        default=False,
        help="""Copy T2-weighted image instead of registering
                    it to the T1-weighted image""",
    ),
)
no_t1 = CommandLineArgument(
    ["--no-t1"],
    dict(
        action="store_true",
        default=False,
        help="""Instructs CHARM that the first image is not a
                    T1-weighted image. In this case, a synthetic T1-like
                    image will be generated from the segmentation and this
                    is used for cortical surface placement (as this
                    operation expects a T1-weighted image).""",
    ),
)
fs_subjects_dir = CommandLineArgument(
    ["--fs-dir"],
    dict(
        type = str,
        help="""If you have a FreeSurfer run of your subject, you can pass the
                    subject directory and CHARM will simply grab the surfaces
                    and spherical registrations from there instead of
                    estimating them. The segmentation will be updated based on
                    these surfaces according to the settings in the .ini file."""
    ),
)
use_settings = CommandLineArgument(
    ["--usesettings"],
    dict(
        nargs=1,
        metavar="settings.ini",
        help="""ini-file with settings (default: charm.ini in
                    simnibs folder)""",
    ),
)
no_neck = CommandLineArgument(
    ["--noneck"],
    dict(
        action="store_true",
        default=False,
        help="""Inform the segmentation that there is no neck in the scan.""",
    ),
)
init_transform = CommandLineArgument(
    ["--inittransform"],
    dict(
        help="""Transformation matrix used
                    to initialize the affine registration of the MNI
                    template to the subject MRI, i.e., it takes the MNI
                    template *to* subject space. Supplied as a path to a
                    space delimited .txt file containing a 4x4
                    transformation matrix (default = None)."""
    ),
)
force_qform = CommandLineArgument(
    ["--forceqform"],
    dict(action="store_true", default=False, help="""Replace sform with qform."""),
)
force_sform = CommandLineArgument(
    ["--forcesform"],
    dict(
        action="store_true",
        default=False,
        help="""Replace qform with sform. Note: strips shears.""",
    ),
)
use_transform = CommandLineArgument(
    ["--usetransform"],
    dict(
        help="""Transformation matrix used
                    instead of doing affine registration of the MNI
                    template to the subject MRI, i.e., it takes the MNI
                    template *to* subject space. Supplied as a path to a
                    space delimited .txt file containing a 4x4
                    transformation matrix (default = None)."""
    ),
)
