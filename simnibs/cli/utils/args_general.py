from simnibs import __version__
from .helpers import CommandLineArgument

subid = CommandLineArgument(
    ["subid"],
    dict(
        type=str,
        help="""Subject ID or /path/to/{subid} or /path/to/m2m_{subid}.
        The former will resolve to "m2m_{subid}" in the current working
        directory. The latter cases both resolve to /path/to/m2m_{subid}.
        """,
    ),
)

version = CommandLineArgument(
    ["-v", "--version"], dict(action="version", version=__version__)
)

debug = CommandLineArgument(
    ['--debug'],
    dict(action='store_true', default=False,
        help="""Write results from intermediate steps to disk."""))