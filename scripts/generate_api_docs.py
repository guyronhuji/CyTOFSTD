"""Generate API reference docs from cytofstandard source code.

This script parses Python source with ``ast`` and writes:

- docs/api/README.md
- docs/api/api_manifest.json
- docs/api/reference/*.md (one file per module)

Run from repository root:

    python3 scripts/generate_api_docs.py
"""

from __future__ import annotations

import ast
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = REPO_ROOT / "cytofstandard"
DOCS_API_DIR = REPO_ROOT / "docs" / "api"
REFERENCE_DIR = DOCS_API_DIR / "reference"
INDEX_PATH = DOCS_API_DIR / "README.md"
MANIFEST_PATH = DOCS_API_DIR / "api_manifest.json"


@dataclass
class FunctionInfo:
    name: str
    signature: str
    docstring: str
    decorators: list[str] = field(default_factory=list)


@dataclass
class ClassInfo:
    name: str
    docstring: str
    bases: list[str] = field(default_factory=list)
    methods: list[FunctionInfo] = field(default_factory=list)


@dataclass
class ModuleInfo:
    module: str
    file_path: str
    docstring: str
    exports: list[str] = field(default_factory=list)
    functions: list[FunctionInfo] = field(default_factory=list)
    classes: list[ClassInfo] = field(default_factory=list)


def short_doc(text: str) -> str:
    if not text:
        return ""
    return text.strip().splitlines()[0].strip()


def annotation_to_text(node: ast.AST | None) -> str:
    if node is None:
        return ""
    try:
        return ast.unparse(node)
    except Exception:
        return ""


def format_arg(arg: ast.arg, default: ast.AST | None = None) -> str:
    name = arg.arg
    annotation = annotation_to_text(arg.annotation)
    if annotation:
        name = f"{name}: {annotation}"
    if default is not None:
        try:
            default_text = ast.unparse(default)
        except Exception:
            default_text = "..."
        name = f"{name} = {default_text}"
    return name


def function_signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    parts: list[str] = []
    args = node.args

    positional = args.posonlyargs + args.args
    defaults = [None] * (len(positional) - len(args.defaults)) + list(args.defaults)

    for arg, default in zip(positional, defaults):
        parts.append(format_arg(arg, default))

    if args.vararg:
        vararg = args.vararg.arg
        vararg_ann = annotation_to_text(args.vararg.annotation)
        if vararg_ann:
            vararg = f"{vararg}: {vararg_ann}"
        parts.append(f"*{vararg}")
    elif args.kwonlyargs:
        parts.append("*")

    for kw_arg, kw_default in zip(args.kwonlyargs, args.kw_defaults):
        parts.append(format_arg(kw_arg, kw_default))

    if args.kwarg:
        kwarg = args.kwarg.arg
        kwarg_ann = annotation_to_text(args.kwarg.annotation)
        if kwarg_ann:
            kwarg = f"{kwarg}: {kwarg_ann}"
        parts.append(f"**{kwarg}")

    result = f"({', '.join(parts)})"
    return_annotation = annotation_to_text(node.returns)
    if return_annotation:
        result += f" -> {return_annotation}"
    return result


def get_decorators(node: ast.FunctionDef | ast.AsyncFunctionDef) -> list[str]:
    decorators: list[str] = []
    for dec in node.decorator_list:
        try:
            decorators.append(ast.unparse(dec))
        except Exception:
            decorators.append("<decorator>")
    return decorators


def parse_exports(tree: ast.Module) -> list[str]:
    exports: list[str] = []
    for stmt in tree.body:
        if not isinstance(stmt, ast.Assign):
            continue
        for target in stmt.targets:
            if isinstance(target, ast.Name) and target.id == "__all__":
                if isinstance(stmt.value, (ast.List, ast.Tuple)):
                    for elt in stmt.value.elts:
                        if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                            exports.append(elt.value)
    return exports


def parse_function(node: ast.FunctionDef | ast.AsyncFunctionDef) -> FunctionInfo:
    return FunctionInfo(
        name=node.name,
        signature=function_signature(node),
        docstring=(ast.get_docstring(node) or "").strip(),
        decorators=get_decorators(node),
    )


def parse_module(path: Path) -> ModuleInfo:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    module_name = "cytofstandard." + str(path.relative_to(PACKAGE_DIR)).replace(
        "/", "."
    ).replace("\\", ".")
    module_name = module_name.removesuffix(".py")
    if module_name.endswith(".__init__"):
        module_name = module_name[: -len(".__init__")]

    info = ModuleInfo(
        module=module_name,
        file_path=str(path.relative_to(REPO_ROOT)),
        docstring=(ast.get_docstring(tree) or "").strip(),
        exports=parse_exports(tree),
    )

    for stmt in tree.body:
        if isinstance(
            stmt, (ast.FunctionDef, ast.AsyncFunctionDef)
        ) and not stmt.name.startswith("_"):
            info.functions.append(parse_function(stmt))
        elif isinstance(stmt, ast.ClassDef) and not stmt.name.startswith("_"):
            bases: list[str] = []
            for base in stmt.bases:
                try:
                    bases.append(ast.unparse(base))
                except Exception:
                    bases.append("<base>")

            cls = ClassInfo(
                name=stmt.name,
                docstring=(ast.get_docstring(stmt) or "").strip(),
                bases=bases,
            )
            for class_stmt in stmt.body:
                if isinstance(
                    class_stmt, (ast.FunctionDef, ast.AsyncFunctionDef)
                ) and not class_stmt.name.startswith("_"):
                    cls.methods.append(parse_function(class_stmt))
            info.classes.append(cls)

    info.functions.sort(key=lambda f: f.name)
    info.classes.sort(key=lambda c: c.name)
    for cls in info.classes:
        cls.methods.sort(key=lambda m: m.name)
    return info


def module_doc_path(module: str) -> Path:
    safe = module.replace(".", "_") + ".md"
    return REFERENCE_DIR / safe


def render_function(func: FunctionInfo, level: str = "###") -> str:
    lines = [f"{level} `{func.name}{func.signature}`"]
    if func.decorators:
        lines.append("")
        lines.append(f"- Decorators: `{', '.join(func.decorators)}`")
    if func.docstring:
        lines.append("")
        lines.append(func.docstring)
    else:
        lines.append("")
        lines.append("No docstring provided.")
    return "\n".join(lines)


def render_module_markdown(info: ModuleInfo) -> str:
    lines: list[str] = []
    lines.append(f"# `{info.module}`")
    lines.append("")
    lines.append(f"- Source: `{info.file_path}`")
    lines.append("")
    if info.docstring:
        lines.append(info.docstring)
        lines.append("")

    if info.exports:
        lines.append("## Public Exports (`__all__`)")
        lines.append("")
        for item in info.exports:
            lines.append(f"- `{item}`")
        lines.append("")

    lines.append("## Top-level Functions")
    lines.append("")
    if info.functions:
        for func in info.functions:
            lines.append(render_function(func, level="###"))
            lines.append("")
    else:
        lines.append("No public top-level functions.")
        lines.append("")

    lines.append("## Classes")
    lines.append("")
    if info.classes:
        for cls in info.classes:
            lines.append(f"### `{cls.name}`")
            lines.append("")
            if cls.bases:
                lines.append(f"- Inherits: `{', '.join(cls.bases)}`")
                lines.append("")
            if cls.docstring:
                lines.append(cls.docstring)
                lines.append("")
            else:
                lines.append("No docstring provided.")
                lines.append("")

            lines.append("#### Methods")
            lines.append("")
            if cls.methods:
                for method in cls.methods:
                    lines.append(render_function(method, level="#####"))
                    lines.append("")
            else:
                lines.append("No public methods.")
                lines.append("")
    else:
        lines.append("No public classes.")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def render_index(modules: list[ModuleInfo], generated_at: str) -> str:
    lines: list[str] = []
    lines.append("# API Reference")
    lines.append("")
    lines.append("Auto-generated from source code. Do not edit these files by hand.")
    lines.append("")
    lines.append(f"- Generated at: `{generated_at}`")
    lines.append(f"- Generator: `scripts/generate_api_docs.py`")
    lines.append("")
    lines.append("## Modules")
    lines.append("")
    lines.append("| Module | Summary | Classes | Functions |")
    lines.append("|---|---|---:|---:|")
    for mod in modules:
        rel_link = module_doc_path(mod.module).relative_to(DOCS_API_DIR)
        summary = short_doc(mod.docstring).replace("|", "\\|")
        lines.append(
            f"| [`{mod.module}`]({rel_link.as_posix()}) | {summary} | {len(mod.classes)} | {len(mod.functions)} |"
        )
    lines.append("")
    lines.append("## Updating")
    lines.append("")
    lines.append("Re-generate docs after any API change:")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/generate_api_docs.py")
    lines.append("```")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_manifest(modules: list[ModuleInfo], generated_at: str) -> dict[str, Any]:
    data: dict[str, Any] = {
        "generated_at": generated_at,
        "generator": "scripts/generate_api_docs.py",
        "package": "cytofstandard",
        "module_count": len(modules),
        "modules": [],
    }

    for mod in modules:
        module_data: dict[str, Any] = {
            "module": mod.module,
            "file_path": mod.file_path,
            "summary": short_doc(mod.docstring),
            "exports": mod.exports,
            "functions": [
                {
                    "name": func.name,
                    "signature": func.signature,
                    "summary": short_doc(func.docstring),
                    "decorators": func.decorators,
                }
                for func in mod.functions
            ],
            "classes": [
                {
                    "name": cls.name,
                    "summary": short_doc(cls.docstring),
                    "bases": cls.bases,
                    "methods": [
                        {
                            "name": method.name,
                            "signature": method.signature,
                            "summary": short_doc(method.docstring),
                            "decorators": method.decorators,
                        }
                        for method in cls.methods
                    ],
                }
                for cls in mod.classes
            ],
        }
        data["modules"].append(module_data)

    return data


def discover_modules() -> list[Path]:
    modules = sorted(
        p for p in PACKAGE_DIR.rglob("*.py") if "__pycache__" not in p.parts
    )
    return modules


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    module_paths = discover_modules()
    modules = [parse_module(path) for path in module_paths]
    modules.sort(key=lambda m: m.module)

    generated_at = datetime.now(timezone.utc).isoformat()

    # Module-level markdown files.
    for mod in modules:
        output_path = module_doc_path(mod.module)
        output_path.write_text(render_module_markdown(mod), encoding="utf-8")

    # Index + manifest.
    INDEX_PATH.write_text(render_index(modules, generated_at), encoding="utf-8")
    MANIFEST_PATH.write_text(
        json.dumps(build_manifest(modules, generated_at), indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Generated API docs for {len(modules)} modules.")
    print(f"- {INDEX_PATH.relative_to(REPO_ROOT)}")
    print(f"- {MANIFEST_PATH.relative_to(REPO_ROOT)}")
    print(f"- {REFERENCE_DIR.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
