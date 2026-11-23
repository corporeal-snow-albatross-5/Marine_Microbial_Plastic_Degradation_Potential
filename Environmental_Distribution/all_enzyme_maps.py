import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
from matplotlib.transforms import offset_copy
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# =========================================================
# 9 ENZYME DATASETS
# =========================================================

DATASETS = [
    {
        "name": "Oxidoreductase (UOM43036.1)",
        "prefix": "oxidoreductase",
        "geo_file": "Oxidoreductase_geographic_locations.csv",
        "meta_file": "Oxidoreductase_isolation_metadata_standardized.xlsx",
        "color": "#FFFB00",  # yellow
    },
    {
        "name": "Poly(3-hydroxybutyrate) depolymerase (LC127088.1)",
        "prefix": "poly3hb_depolymerase",
        "geo_file": "Poly(3-hydroxybutyrate)Depolymerase_geographic_locations.csv",
        "meta_file": "Poly(3-hydroxybutyrate)Depolymerase_isolation_metadata_standardized.xlsx",
        "color": "#009E3F",  # green
    },
    {
        "name": "Glutathione peroxidase (WP_026826575.1)",
        "prefix": "gpx_wp_026826575_1",
        "geo_file": "GlutathionePeroxidase(WP_026826575.1)_geographic_locations.csv",
        "meta_file": "GlutathionePeroxidase(WP_026826575.1)_isolation_metadata_standardized.xlsx",
        "color": "#0072B2",  # blue
    },
    {
        "name": "Glutathione peroxidase (WP_003945816.1)",
        "prefix": "gpx_wp_003945816_1",
        "geo_file": "GlutathionePeroxidase(WP_003945816.1)_geographic_locations.csv",
        "meta_file": "GlutathionePeroxidase(WP_003945816.1)_isolation_metadata_standardized.xlsx",
        "color": "#42d4f4",  # cyan
    },
    {
        "name": "PHB depolymerase (AAB40611.1)",
        "prefix": "phb",
        "geo_file": "PHB_geographic_locations.csv",
        "meta_file": "PHB_isolation_metadata_standardized.xlsx",
        "color": "#8E63CE",  # purple
    },
    {
        "name": "Laccase (UVG67878.1)",
        "prefix": "laccase",
        "geo_file": "Laccase_geographic_locations.csv",
        "meta_file": "Laccase_isolation_metadata_standardized.xlsx",
        "color": "#7B3E1B",  # brown
    },
    {
        "name": "MCL PHA depolymerase (AAQ72538.1)",
        "prefix": "pha",
        "geo_file": "PHA_geographic_locations.csv",
        "meta_file": "PHA_isolation_metadata_standardized.xlsx",
        "color": "#F079BA",  # pink
    },
    {
        "name": "PETase (UUT36763.1)",
        "prefix": "petase184",
        "geo_file": "PETase184_geographic_locations.csv",
        "meta_file": "PETase184_metadata_standardized.xlsx",
        "color": "#FF0808",  # red
    },
    {
        "name": "PETase (UUT36764.1)",
        "prefix": "petase183",
        "geo_file": "PETase183_geographic_locations.csv",
        "meta_file": "PETase183_metadata_standardized.xlsx",
        "color": "#E69F00",  # orange
    },
]


# =========================================================
# ENVIRONMENTAL GROUP SYMBOLS
# =========================================================

ENV_STYLES = {
    "Plastisphere / marine biofilm✕": {
        "label": "Plastisphere / marine biofilm",
        "marker": "X",
    },

    "Marine host-associated⬡": {
        "label": "Marine host-associated",
        "marker": "h",
    },

    "Mangrove sediment✚": {
        "label": "Mangrove sediment",
        "marker": "P",
    },

    "Polar sea ice / cryosphere◆": {
        "label": "Polar sea ice / cryosphere",
        "marker": "D",
    },

    "Deep-sea / chemosynthetic sediment⬟": {
        "label": "Deep-sea / chemosynthetic sediment",
        "marker": "p",
    },

    "Intertidal sediment / beach sand▲": {
        "label": "Intertidal sediment / beach sand",
        "marker": "^",
    },

    "Intertidal / coastal seawater★": {
        "label": "Intertidal / coastal seawater",
        "marker": "*",
    },

    "Deep seawater▼": {
        "label": "Deep seawater",
        "marker": "v",
    },

    "Marine sediment■": {
        "label": "Marine sediment",
        "marker": "s",
    },

    "Surface seawater●": {
        "label": "Surface seawater",
        "marker": "o",
    },
}


ENV_ORDER = [
    "Plastisphere / marine biofilm✕",
    "Marine host-associated⬡",
    "Mangrove sediment✚",
    "Polar sea ice / cryosphere◆",
    "Deep-sea / chemosynthetic sediment⬟",
    "Intertidal sediment / beach sand▲",
    "Intertidal / coastal seawater★",
    "Deep seawater▼",
    "Marine sediment■",
    "Surface seawater●",
]


# =========================================================
# MARKER SCALE
# =========================================================

MARKER_SCALE = {
    "X": 0.85,
    "h": 1.00,
    "P": 0.90,
    "D": 0.85,
    "p": 1.00,
    "^": 1.00,
    "*": 1.00,
    "v": 1.00,
    "s": 0.95,
    "o": 1.00
}


BASE_SIZE = 6

# Repeated samples still get larger, but growth is deliberately
# sub-linear so highly repeated sites do not become enormous.
REPEAT_INCREMENT = 2.0

# Jitter is applied only when multiple distinct plotted markers share
# essentially the same geographic coordinate. Offsets are in display
# points, so they are visual only and do not alter the stored lat/lon.
JITTER_COORD_DECIMALS = 3
JITTER_BASE_RADIUS_PT = 5.0
JITTER_RADIUS_STEP_PT = 0.65
JITTER_MAX_RADIUS_PT = 11.0

# Transparent enough to see overlap while retaining clear color
FILL_ALPHA = 0.45
EDGE_ALPHA = 0.90


# =========================================================
# MAP STYLE
# =========================================================

OCEAN_COLOR = "#cfe3f1"
LAND_COLOR = "#ffffff"

BORDER_COLOR = "#aaaaaa"
COAST_COLOR = "#969696"

LAND_LABEL_COLOR = "#444444"
OCEAN_LABEL_COLOR = "#3168b3"

LABEL_SIZE = 12


# =========================================================
# LABEL LOCATIONS
# =========================================================

CONTINENT_LABELS = [
    ("North America", -95, 48),
    ("South America", -58, -8),
    ("Europe", 17, 50),
    ("Africa", 20, 5),
    ("Asia", 90, 46),
    ("Australia", 134, -26),
    ("Antarctica", 25, -78)
]


OCEAN_LABELS = [
    ("Pacific Ocean", -145, 8),
    ("Pacific Ocean", 155, 8),
    ("Atlantic Ocean", -45, 19),
    ("Indian Ocean", 77, -24),
    ("Arctic Ocean", 3, 73),
    ("Southern Ocean", 25, -57)
]


# =========================================================
# FALLBACK COORDINATES
#
# Only used if no numeric coordinates can be found.
# =========================================================

LOCATION_HINTS = {
    "atlantic ocean at key west": (24.55, -81.78),
    "florida keys": (24.55, -81.78),

    "bagamoyo": (-6.44, 38.90),

    "matang mangrove perak": (4.84, 100.63),

    "algoa bay": (-33.96, 25.64),
    "port elizabeth": (-33.96, 25.64),

    "arctic ocean": (79.80, 1.76),

    "south china sea": (18.00, 116.00),
    "nanhai": (18.00, 116.00),

    "taiwan": (23.70, 121.00),

    "black sea": (43.00, 34.00),

    "florianopolis": (-27.60, -48.55),

    "mariana trench": (11.40, 142.40),

    "weizhou island": (21.03, 109.12),

    "beigang island": (21.58, 109.10),

    "quang ngai": (15.12, 108.80),

    "sea of japan": (43.50, 135.50),

    "adelaide island": (-67.34, -68.08),

    "xisha": (16.75, 112.35),
}


# =========================================================
# READ METADATA
# =========================================================

def read_metadata(meta_file):

    meta = pd.read_excel(
        meta_file,
        sheet_name=0
    )

    keep = [
        "Protein accession number",
        "Environment group"
    ]

    if "Isolation location" in meta.columns:
        keep.append("Isolation location")

    return meta[keep].copy()


# =========================================================
# SPLIT CELLS CONTAINING MULTIPLE VALUES
# =========================================================

def split_multivalue(value):

    if pd.isna(value):
        return []

    if isinstance(
        value,
        (int, float, np.integer, np.floating)
    ):
        return [value]

    text = str(value).strip()

    if not text:
        return []

    if text.lower().startswith("not reported"):
        return []

    pieces = [
        piece.strip()
        for piece in text.split("|")
    ]

    cleaned = []

    for piece in pieces:

        piece = re.sub(
            r"^[A-Z]{2,4}_[^:]+:\s*",
            "",
            piece
        )

        cleaned.append(
            piece.strip()
        )

    return cleaned


# =========================================================
# CONVERT ONE COORDINATE VALUE TO NUMBER
# =========================================================

def parse_coordinate_token(value):

    if pd.isna(value):
        return np.nan

    if isinstance(
        value,
        (int, float, np.integer, np.floating)
    ):
        return float(value)

    text = str(value).strip()

    if not text:
        return np.nan

    if text.lower().startswith("not reported"):
        return np.nan

    match = re.search(
        r"(-?\d+(?:\.\d+)?)",
        text
    )

    if not match:
        return np.nan

    number = float(
        match.group(1)
    )

    direction = re.search(
        r"([NSEW])\b",
        text,
        re.I
    )

    if direction:

        if direction.group(1).upper() in (
            "S",
            "W"
        ):
            number = -abs(number)

    return number


# =========================================================
# FIND COORDINATES EMBEDDED IN TEXT
# =========================================================

def find_coordinate_pairs_in_text(text):

    if pd.isna(text):
        return []

    pattern = re.compile(
        r"(-?\d+(?:\.\d+)?)\s*°?\s*([NS])\s*[,;\s]+"
        r"(-?\d+(?:\.\d+)?)\s*°?\s*([EW])",
        re.I,
    )

    pairs = []

    for match in pattern.finditer(
        str(text)
    ):

        lat = float(
            match.group(1)
        )

        lon = float(
            match.group(3)
        )

        if match.group(2).upper() == "S":
            lat = -abs(lat)

        if match.group(4).upper() == "W":
            lon = -abs(lon)

        pairs.append(
            (lat, lon)
        )

    return pairs


# =========================================================
# APPROXIMATE A NAMED LOCATION
# =========================================================

def approximate_from_location(text):

    if pd.isna(text):
        return None

    text = str(text).lower()

    for location_name, coordinates in LOCATION_HINTS.items():

        if location_name in text:

            return coordinates

    return None


# =========================================================
# EXPAND EACH DATASET INTO PLOTTABLE POINTS
# =========================================================

def expand_points(
    merged_df,
    enzyme_name,
    color
):

    mapped = []
    unmapped = []

    text_fields = [
        "Exact sampling location",
        "Coordinate evidence",
        "Isolation location",
        "Location evidence",
    ]

    for _, row in merged_df.iterrows():

        points = []

        # ---------------------------------------------
        # 1. Numeric latitude / longitude columns
        # ---------------------------------------------

        lat_parts = split_multivalue(
            row.get("Latitude")
        )

        lon_parts = split_multivalue(
            row.get("Longitude")
        )

        if (
            lat_parts
            and lon_parts
            and len(lat_parts) == len(lon_parts)
        ):

            for lat_part, lon_part in zip(
                lat_parts,
                lon_parts
            ):

                lat = parse_coordinate_token(
                    lat_part
                )

                lon = parse_coordinate_token(
                    lon_part
                )

                if (
                    not np.isnan(lat)
                    and not np.isnan(lon)
                ):

                    points.append(
                        (lat, lon)
                    )

        else:

            lat = parse_coordinate_token(
                row.get("Latitude")
            )

            lon = parse_coordinate_token(
                row.get("Longitude")
            )

            if (
                not np.isnan(lat)
                and not np.isnan(lon)
            ):

                points.append(
                    (lat, lon)
                )


        # ---------------------------------------------
        # 2. Coordinates written inside text
        # ---------------------------------------------

        if not points:

            for field in text_fields:

                points.extend(
                    find_coordinate_pairs_in_text(
                        row.get(field)
                    )
                )


        # ---------------------------------------------
        # 3. Approximate named sampling location
        # ---------------------------------------------

        if not points:

            for field in text_fields:

                approximate = (
                    approximate_from_location(
                        row.get(field)
                    )
                )

                if approximate is not None:

                    points = [
                        approximate
                    ]

                    break


        # ---------------------------------------------
        # Remove duplicate coordinates within record
        # ---------------------------------------------

        unique_points = []
        seen = set()

        for lat, lon in points:

            key = (
                round(float(lat), 5),
                round(float(lon), 5)
            )

            if key not in seen:

                unique_points.append(
                    (
                        float(lat),
                        float(lon)
                    )
                )

                seen.add(key)


        record = {

            "enzyme":
                enzyme_name,

            "color":
                color,

            "Protein accession number":
                row.get(
                    "Protein accession number"
                ),

            "Environment group":
                row.get(
                    "Environment group"
                ),

            "Exact sampling location":
                row.get(
                    "Exact sampling location",
                    row.get(
                        "Isolation location",
                        ""
                    )
                ),
        }


        if unique_points:

            for lat, lon in unique_points:

                mapped.append(
                    {
                        **record,
                        "lat": lat,
                        "lon": lon
                    }
                )

        else:

            unmapped.append(
                record
            )


    return (
        pd.DataFrame(mapped),
        pd.DataFrame(unmapped)
    )


# =========================================================
# CREATE BASE MAP
# =========================================================

def add_base_map(ax):

    ax.set_global()

    ax.add_feature(
        cfeature.OCEAN,
        facecolor=OCEAN_COLOR,
        zorder=0
    )

    ax.add_feature(
        cfeature.LAND,
        facecolor=LAND_COLOR,
        zorder=1
    )

    ax.add_feature(
        cfeature.BORDERS,
        edgecolor=BORDER_COLOR,
        linewidth=0.50,
        zorder=2
    )

    ax.coastlines(
        color=COAST_COLOR,
        linewidth=0.65,
        zorder=2
    )


    # ---------------------------------------------
    # Continents
    # ---------------------------------------------

    for text, lon, lat in CONTINENT_LABELS:

        ax.text(
            lon,
            lat,
            text,

            transform=ccrs.PlateCarree(),

            ha="center",
            va="center",

            fontsize=LABEL_SIZE,

            color=LAND_LABEL_COLOR,

            fontweight="bold",

            zorder=3,

            clip_on=True
        )


    # ---------------------------------------------
    # Oceans
    # ---------------------------------------------

    for text, lon, lat in OCEAN_LABELS:

        ax.text(
            lon,
            lat,
            text,

            transform=ccrs.PlateCarree(),

            ha="center",
            va="center",

            fontsize=LABEL_SIZE,

            color=OCEAN_LABEL_COLOR,

            fontstyle="italic",

            zorder=3,

            clip_on=True
        )


# =========================================================
# ENVIRONMENT SYMBOL LEGEND
# =========================================================

def environment_legend_handles():

    face = to_rgba(
        "#666666",
        0.55
    )

    edge = to_rgba(
        "#555555",
        0.95
    )

    handles = []

    for env_group in ENV_ORDER:

        marker = (
            ENV_STYLES[env_group]["marker"]
        )

        label = (
            ENV_STYLES[env_group]["label"]
        )

        size = (
            9
            * MARKER_SCALE.get(
                marker,
                1.0
            )
        )

        handles.append(

            Line2D(
                [0],
                [0],

                marker=marker,
                linestyle="None",

                markersize=size,

                markerfacecolor=face,
                markeredgecolor=edge,

                markeredgewidth=1.15,

                label=label
            )
        )

    return handles


# =========================================================
# ENZYME COLOR LEGEND
# =========================================================

def enzyme_legend_handles():

    handles = []

    for config in DATASETS:

        face = to_rgba(
            config["color"],
            FILL_ALPHA
        )

        edge = to_rgba(
            config["color"],
            EDGE_ALPHA
        )

        handles.append(

            Line2D(
                [0],
                [0],

                marker="o",
                linestyle="None",

                markersize=8,

                markerfacecolor=face,
                markeredgecolor=edge,

                markeredgewidth=1.15,

                label=config["name"]
            )
        )

    return handles


# =========================================================
# REPEAT-SIZE + JITTER HELPERS
# =========================================================

def marker_size_for_repeat(repeat_count, marker):

    # Diminishing growth: every additional repeat still increases
    # the symbol, but by less than the previous linear +4 rule.
    repeats_above_one = max(float(repeat_count) - 1.0, 0.0)

    raw_size = (
        BASE_SIZE
        + REPEAT_INCREMENT * np.sqrt(repeats_above_one)
    )

    return (
        raw_size
        * MARKER_SCALE.get(marker, 1.0)
    )


def add_jitter_offsets(grouped_df):

    jittered = grouped_df.copy()

    jittered["jitter_x_pt"] = 0.0
    jittered["jitter_y_pt"] = 0.0

    # Rounding makes coordinates that differ only at tiny numerical
    # precision count as the same plotted location.
    jittered["_jitter_lat_key"] = (
        jittered["lat"].astype(float).round(JITTER_COORD_DECIMALS)
    )

    jittered["_jitter_lon_key"] = (
        jittered["lon"].astype(float).round(JITTER_COORD_DECIMALS)
    )

    grouped_indices = jittered.groupby(
        ["_jitter_lat_key", "_jitter_lon_key"],
        sort=False
    ).groups

    for _, indices in grouped_indices.items():

        indices = list(indices)
        n = len(indices)

        if n <= 1:
            continue

        # Stable order means the same dataset always gets the same
        # jitter arrangement every time the script is run.
        ordered_indices = (
            jittered.loc[indices]
            .sort_values(
                ["enzyme", "Environment group"],
                kind="stable"
            )
            .index
            .tolist()
        )

        radius = min(
            JITTER_BASE_RADIUS_PT
            + JITTER_RADIUS_STEP_PT * (n - 2),
            JITTER_MAX_RADIUS_PT
        )

        for position, row_index in enumerate(ordered_indices):

            angle = (
                np.pi / 2
                + 2 * np.pi * position / n
            )

            jittered.at[row_index, "jitter_x_pt"] = (
                radius * np.cos(angle)
            )

            jittered.at[row_index, "jitter_y_pt"] = (
                radius * np.sin(angle)
            )

    return jittered.drop(
        columns=["_jitter_lat_key", "_jitter_lon_key"]
    )


# =========================================================
# PLOT ONE MAP
# =========================================================

def plot_map(
    grouped_df,
    title,
    png_path,
    pdf_path,
    combined=False
):

    fig = plt.figure(
        figsize=(16, 10)
    )

    ax = plt.axes(
        projection=ccrs.Robinson()
    )

    add_base_map(ax)


    fig.subplots_adjust(
        bottom=0.34 if combined else 0.23,
        top=0.92
    )


    # ---------------------------------------------
    # Plot sampling points
    # ---------------------------------------------

    # Distinct markers at the same coordinate are separated slightly
    # in display space. Isolated markers remain exactly where they are.
    plot_df = add_jitter_offsets(grouped_df)

    plate_carree_transform = (
        ccrs.PlateCarree()._as_mpl_transform(ax)
    )

    for _, row in plot_df.iterrows():

        env_group = (
            row["Environment group"]
        )

        marker = ENV_STYLES.get(
            env_group,
            {"marker": "o"}
        )["marker"]


        marker_size = marker_size_for_repeat(
            row["repeat_count"],
            marker
        )


        face = to_rgba(
            row["color"],
            FILL_ALPHA
        )

        edge = to_rgba(
            row["color"],
            EDGE_ALPHA
        )


        jitter_x = float(row["jitter_x_pt"])
        jitter_y = float(row["jitter_y_pt"])

        if jitter_x != 0.0 or jitter_y != 0.0:

            point_transform = offset_copy(
                plate_carree_transform,
                fig=fig,
                x=jitter_x,
                y=jitter_y,
                units="points"
            )

        else:

            point_transform = ccrs.PlateCarree()


        ax.plot(
            row["lon"],
            row["lat"],

            marker=marker,
            linestyle="None",

            markersize=marker_size,

            markerfacecolor=face,
            markeredgecolor=edge,

            markeredgewidth=1.15,

            transform=point_transform,

            zorder=10
        )


    # ---------------------------------------------
    # Title
    # ---------------------------------------------

    plt.title(
        title,

        fontsize=20,

        fontweight="bold",

        pad=14
    )


    # ---------------------------------------------
    # Environmental-group legend
    # ---------------------------------------------

    env_legend = fig.legend(

        handles=
            environment_legend_handles(),

        loc="lower center",

        bbox_to_anchor=(
            0.5,
            0.03 if not combined else 0.035
        ),

        ncol=
            2 if not combined else 5,

        frameon=True,

        title=
            "Environmental group = symbol",

        fontsize=10,

        title_fontsize=11,

        handletextpad=0.7,

        columnspacing=1.4,

        borderpad=0.8
    )


    env_legend.get_frame().set_facecolor(
        "white"
    )

    env_legend.get_frame().set_edgecolor(
        "#c8c8c8"
    )

    env_legend.get_frame().set_alpha(
        0.97
    )


    # ---------------------------------------------
    # Combined map also gets enzyme color legend
    # ---------------------------------------------

    if combined:

        fig.add_artist(
            env_legend
        )

        enzyme_legend = fig.legend(

            handles=
                enzyme_legend_handles(),

            loc="lower center",

            bbox_to_anchor=(
                0.5,
                0.19
            ),

            ncol=3,

            frameon=True,

            title="Enzyme = color",

            fontsize=10,

            title_fontsize=11,

            handletextpad=0.7,

            columnspacing=1.4,

            borderpad=0.8
        )


        enzyme_legend.get_frame().set_facecolor(
            "white"
        )

        enzyme_legend.get_frame().set_edgecolor(
            "#c8c8c8"
        )

        enzyme_legend.get_frame().set_alpha(
            0.97
        )


    # ---------------------------------------------
    # Save
    # ---------------------------------------------

    plt.savefig(
        png_path,

        dpi=400,

        bbox_inches="tight"
    )


    plt.savefig(
        pdf_path,

        bbox_inches="tight"
    )


    plt.close(fig)


# =========================================================
# MAIN
# =========================================================

def main():

    base_dir = os.getcwd()

    output_dir = os.path.join(
        base_dir,
        "outputs"
    )

    individual_dir = os.path.join(
        output_dir,
        "individual_maps"
    )

    os.makedirs(
        individual_dir,
        exist_ok=True
    )


    all_grouped = []
    all_mapped = []
    all_unmapped = []


    # =====================================================
    # PROCESS EACH ENZYME
    # =====================================================

    for config in DATASETS:

        geo_path = os.path.join(
            base_dir,
            config["geo_file"]
        )

        meta_path = os.path.join(
            base_dir,
            config["meta_file"]
        )


        # Check files exist
        if not os.path.exists(geo_path):

            raise FileNotFoundError(
                f"Missing file: "
                f"{config['geo_file']}"
            )


        if not os.path.exists(meta_path):

            raise FileNotFoundError(
                f"Missing file: "
                f"{config['meta_file']}"
            )


        # ---------------------------------------------
        # Read
        # ---------------------------------------------

        geo = pd.read_csv(
            geo_path
        )

        meta = read_metadata(
            meta_path
        )


        # ---------------------------------------------
        # Merge geographic + metadata tables
        # ---------------------------------------------

        merged = geo.merge(

            meta,

            on=
                "Protein accession number",

            how="left"
        )


        # ---------------------------------------------
        # Convert rows into sampling points
        # ---------------------------------------------

        mapped, unmapped = expand_points(

            merged,

            config["name"],

            config["color"]
        )


        # ---------------------------------------------
        # Group repeated sampling locations
        #
        # Same:
        # enzyme
        # + location
        # + environmental group
        #
        # becomes one larger marker.
        # ---------------------------------------------

        if not mapped.empty:

            grouped = (

                mapped.groupby(

                    [
                        "enzyme",
                        "color",
                        "lat",
                        "lon",
                        "Environment group"
                    ],

                    as_index=False
                )

                .size()

                .rename(
                    columns={
                        "size":
                            "repeat_count"
                    }
                )
            )


            # -----------------------------------------
            # Individual enzyme map
            # -----------------------------------------

            plot_map(

                grouped,

                (
                    "Global Distribution of "
                    f"{config['name']} "
                    "Isolation Sites"
                ),

                os.path.join(
                    individual_dir,
                    (
                        f"{config['prefix']}"
                        "_global_map.png"
                    )
                ),

                os.path.join(
                    individual_dir,
                    (
                        f"{config['prefix']}"
                        "_global_map.pdf"
                    )
                ),

                combined=False
            )


            all_grouped.append(
                grouped
            )

            all_mapped.append(
                mapped
            )


        if not unmapped.empty:

            all_unmapped.append(
                unmapped
            )


        print(
            f"{config['name']}: "
            f"mapped {len(mapped)} point rows, "
            f"unmapped {len(unmapped)} rows"
        )


    # =====================================================
    # SAVE AUDIT DATA
    # =====================================================

    if all_mapped:

        pd.concat(
            all_mapped,
            ignore_index=True
        ).to_csv(

            os.path.join(
                output_dir,
                "all_mapped_points.csv"
            ),

            index=False
        )


    if all_unmapped:

        pd.concat(
            all_unmapped,
            ignore_index=True
        ).to_csv(

            os.path.join(
                output_dir,
                "all_unmapped_rows.csv"
            ),

            index=False
        )

    else:

        pd.DataFrame().to_csv(

            os.path.join(
                output_dir,
                "all_unmapped_rows.csv"
            ),

            index=False
        )


    # =====================================================
    # COMBINED MAP
    # =====================================================

    combined = pd.concat(
        all_grouped,
        ignore_index=True
    )


    combined.to_csv(

        os.path.join(
            output_dir,
            "combined_grouped_points.csv"
        ),

        index=False
    )


    plot_map(

        combined,

        (
            "Global Distribution of Marine "
            "Plastic-Degrading Enzyme Isolation Sites"
        ),

        os.path.join(
            output_dir,
            "combined_enzyme_global_map.png"
        ),

        os.path.join(
            output_dir,
            "combined_enzyme_global_map.pdf"
        ),

        combined=True
    )


    # =====================================================
    # DONE
    # =====================================================

    print()

    print("DONE")

    print(
        "Outputs saved in ./outputs"
    )

    print(
        "  individual_maps/"
        "       9 individual maps"
    )

    print(
        "  combined_enzyme_global_map.png"
    )

    print(
        "  combined_enzyme_global_map.pdf"
    )

    print(
        "  all_mapped_points.csv"
    )

    print(
        "  all_unmapped_rows.csv"
    )

    print(
        "  combined_grouped_points.csv"
    )


if __name__ == "__main__":
    main()