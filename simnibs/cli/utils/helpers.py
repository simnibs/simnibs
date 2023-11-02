from dataclasses import dataclass
from pathlib import Path


def resolve_subject_id_path(subid):
    m2m_dir = Path(subid).resolve()
    return (
        m2m_dir
        if m2m_dir.name.startswith("m2m_")
        else m2m_dir.with_name(f"m2m_{m2m_dir.name}")
    )


# Arguments

@dataclass()
class CommandLineArgument:
    flags: list  # [str]
    actions: dict

    def add_to(self, parser):
        parser.add_argument(*self.flags, **self.actions)


def add_argument(parser, argument: CommandLineArgument):
    parser.add_argument(*argument.flags, **argument.actions)


# Parser

@dataclass()
class CommandLineParser:
    name: str
    kwargs: dict

    def add_to(self, parser, parents = None):
        if parents is None:
            parents = []
        return parser.add_parser(self.name, parents=parents, **self.kwargs)
