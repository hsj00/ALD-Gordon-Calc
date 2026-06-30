"""export.py — 결과/그래프 내보내기 헬퍼 (CSV, PNG)."""
import io
import csv
import datetime


def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _g(x):
    """숫자는 6 유효숫자, 그 외 문자열."""
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, (int, float)):
        return f"{x:.6g}"
    return str(x)


def summary_csv(rows):
    """rows: list of (parameter, value, unit). Excel 한글 호환 BOM 포함."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["parameter", "value", "unit"])
    for p, v, unit in rows:
        w.writerow([p, _g(v), unit])
    return buf.getvalue().encode("utf-8-sig")


def series_csv(colnames, columns):
    """colnames: list[str], columns: list of equal-length sequences."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(colnames)
    for row in zip(*columns):
        w.writerow([_g(v) for v in row])
    return buf.getvalue().encode("utf-8-sig")


def fig_png_bytes(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    return buf.getvalue()
