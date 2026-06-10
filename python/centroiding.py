"""
====================================================
CENTROIDING — Centroid Detection Algorithms
====================================================

Purpose:
    Detects the (x, y) pixel coordinates of stars in a grayscale image.
    Also provides image thresholding methods used to separate stars
    from background noise.

    Contains:
        - BasicThreshold, OtsusThreshold (image thresholding)
        - CogCentroider (Center of Gravity centroiding)
        - IterativeWeightedCogCentroider (Iterative Weighted CoG)
        - DummyCentroider (debug centroids)

    All centroiders follow the same interface:
        Go(image, width, height) -> list[Star]

Thresholding Background:
    Star images have a dark background with bright star spots.
    A threshold separates "likely star pixels" from "likely background".
    BasicThreshold: mean + 5*stddev (works well for star images)
    Otsu: maximizes inter-class variance (general-purpose)

CoG Background:
    The Center of Gravity method finds centroid by intensity-weighted
    average of pixel coordinates within a connected bright region:
        x_c = Σ(x_i * I_i) / Σ(I_i)
        y_c = Σ(y_i * I_i) / Σ(I_i)

    Regions touching the image edge are rejected (they are partial stars).

Reference:
    Original C++ source: src/centroiders.hpp, src/centroiders.cpp
"""

from __future__ import annotations

import math
import random
from typing import List, Set

from math_utils import Vec2
from database import Star, Stars


# ─────────────────────────────────────────
# THRESHOLDING ALGORITHMS
# ─────────────────────────────────────────


def basic_threshold(image: bytes, image_width: int, image_height: int) -> int:
    total_pixels = image_height * image_width
    total_mag = 0
    for i in range(total_pixels):
        total_mag += image[i]

    # C++ uses integer division: totalMag / totalPixels (both integral types)
    mean = total_mag // total_pixels

    std_sum = 0.0
    for i in range(total_pixels):
        diff = image[i] - mean
        std_sum += diff * diff

    std = math.sqrt(std_sum / total_pixels)
    return int(mean + std * 5)


def basic_threshold_onepass(image: bytes, image_width: int,
                            image_height: int) -> int:
    total_pixels = image_height * image_width
    total_mag = 0
    sq_total_mag = 0

    for i in range(total_pixels):
        total_mag += image[i]
        sq_total_mag += image[i] * image[i]

    # C++ uses integer division for the mean
    mean = total_mag / total_pixels
    variance = (sq_total_mag / total_pixels) - (mean * mean)
    std = math.sqrt(variance)
    return int(mean + std * 5)


def otsus_threshold(image: bytes, image_width: int, image_height: int) -> int:
    """Otsu's thresholding method — maximizes inter-class variance.

    Not specifically tailored to star images, but works as a
    general-purpose threshold.

    Reference:
        Original C++ source: src/centroiders.cpp:41-78
    """
    total = image_width * image_height
    histogram = [0] * 256

    for i in range(total):
        histogram[image[i]] += 1

    sum_total = 0.0
    for i in range(256):
        sum_total += i * histogram[i]

    sum_b = 0.0
    w_b = 0.0
    maximum = 0.0
    level = 0

    for i in range(256):
        w_f = total - w_b
        if w_b > 0 and w_f > 0:
            m_f = (sum_total - sum_b) / w_f
            val = w_b * w_f * ((sum_b / w_b) - m_f) * ((sum_b / w_b) - m_f)
            if val >= maximum:
                level = i
                maximum = val
        w_b += histogram[i]
        sum_b += i * histogram[i]

    return level


# ─────────────────────────────────────────
# CENTROIDER BASE CLASS
# ─────────────────────────────────────────


class CentroidAlgorithm:
    """Base class for centroid detection algorithms.

    An algorithm that detects bright points in an image and returns
    their (x, y) coordinates as Star objects.

    Reference:
        Original C++ source: src/centroiders.hpp:12-22
    """

    def go(self, image: bytes, image_width: int, image_height: int) -> Stars:
        """Detect centroids in a grayscale image.

        Args:
            image: Row-major array of grayscale pixels (0-255).
            image_width: Image width in pixels.
            image_height: Image height in pixels.

        Returns:
            List of detected Star objects with pixel coordinates.
        """
        raise NotImplementedError


# ─────────────────────────────────────────
# DUMMY CENTROIDER (debug)
# ─────────────────────────────────────────


class DummyCentroider(CentroidAlgorithm):
    """A centroid algorithm for debugging that returns random centroids.

    Reference:
        Original C++ source: src/centroiders.cpp:19-28
    """

    def __init__(self, num_stars: int = 5):
        self.num_stars = num_stars

    def go(self, image: bytes, image_width: int, image_height: int) -> Stars:
        rng = random.Random(123456)
        result = []
        for _ in range(self.num_stars):
            x = rng.randint(0, image_width - 1)
            y = rng.randint(0, image_height - 1)
            result.append(Star(x, y, 10.0))
        return result


# ─────────────────────────────────────────
# CENTER OF GRAVITY CENTROIDER
# ─────────────────────────────────────────


class CogCentroider(CentroidAlgorithm):
    """Center of Gravity centroid detection.

    Algorithm:
        1. Compute threshold using BasicThreshold.
        2. Scan pixels; for each pixel ≥ threshold not yet visited:
           a. Flood-fill to find all connected bright pixels.
           b. Reject region if it touches the image edge.
           c. Compute centroid as intensity-weighted average of coordinates.
           d. Compute radius from bounding box.
        3. Return list of centroids with +0.5 pixel offset
           (pixel-center convention).

    Complexity: O(width × height × region_area) worst case.

    Reference:
        Original C++ source: src/centroiders.cpp:158-196
    """

    def go(self, image: bytes, image_width: int, image_height: int) -> Stars:
        cutoff = basic_threshold(image, image_width, image_height)
        total_pixels = image_height * image_width
        checked: Set[int] = set()
        result: Stars = []

        for i in range(total_pixels):
            if image[i] < cutoff or i in checked:
                continue

            # Flood-fill to collect all pixels in this star region
            y_coord_sum = 0.0
            x_coord_sum = 0.0
            mag_sum = 0
            x_min = i % image_width
            x_max = x_min
            y_min = i // image_width
            y_max = y_min
            is_valid = True
            size_before = len(checked)

            # Iterative stack-based flood-fill (avoids Python recursion limits)
            stack = [i]
            while stack:
                idx = stack.pop()
                if idx < 0 or idx >= total_pixels:
                    continue
                if image[idx] < cutoff or idx in checked:
                    continue

                # Check if pixel is on the edge of the image
                if (idx % image_width == 0
                        or idx % image_width == image_width - 1
                        or idx // image_width == 0
                        or idx // image_width == image_height - 1):
                    is_valid = False

                checked.add(idx)

                px = idx % image_width
                py = idx // image_width
                if px > x_max:
                    x_max = px
                elif px < x_min:
                    x_min = px
                if py > y_max:
                    y_max = py
                elif py < y_min:
                    y_min = py

                mag_sum += image[idx]
                x_coord_sum += px * image[idx]
                y_coord_sum += py * image[idx]

                # Add neighbors (4-connected) in reverse order so RIGHT is
                # popped first (LIFO), matching C++ recursive right-first traversal
                stack.append(idx - image_width)  # UP
                stack.append(idx + image_width)  # DOWN
                if px != 0:
                    stack.append(idx - 1)        # LEFT
                if px != image_width - 1:
                    stack.append(idx + 1)        # RIGHT — last pushed = first popped

            x_diameter = (x_max - x_min) + 1
            y_diameter = (y_max - y_min) + 1

            if mag_sum > 0 and is_valid:
                x_center = x_coord_sum / mag_sum + 0.5
                y_center = y_coord_sum / mag_sum + 0.5
                result.append(Star(
                    x_center, y_center,
                    x_diameter / 2.0, y_diameter / 2.0,
                    len(checked) - size_before,
                ))

        return result


# ─────────────────────────────────────────
# ITERATIVE WEIGHTED CENTER OF GRAVITY
# ─────────────────────────────────────────


class IterativeWeightedCogCentroider(CentroidAlgorithm):
    """Iterative Weighted Center of Gravity centroid detection.

    A more sophisticated centroid algorithm using Gaussian-weighted
    iterative refinement. Despite the additional complexity, LOST's
    authors report it doesn't perform much better than CoG.

    Algorithm:
        1. Compute threshold via BasicThreshold.
        2. Flood-fill to find connected bright region.
        3. Find pixel with maximum intensity as initial guess.
        4. Compute Full Width at Half Maximum (FWHM).
        5. Convert FWHM to Gaussian standard deviation:
           σ = FWHM / (2 * sqrt(2 * ln(2)))
        6. Iteratively re-estimate centroid using Gaussian weights:
           w = I_max * exp(-((x-x_g)² + (y-y_g)²) / (2σ²))
        7. Stop when centroid change < 0.0002 or 100k iterations.

    Reference:
        Original C++ source: src/centroiders.cpp:247-325
    """

    MIN_CHANGE = 0.0002  # Convergence threshold

    def go(self, image: bytes, image_width: int, image_height: int) -> Stars:
        cutoff = basic_threshold(image, image_width, image_height)
        total_pixels = image_height * image_width
        checked: Set[int] = set()
        result: Stars = []

        for i in range(total_pixels):
            if image[i] < cutoff or i in checked:
                continue

            # Flood-fill to collect star indices
            star_indices: List[int] = []
            max_intensity = 0
            guess = i
            x_min = i % image_width
            x_max = x_min
            y_min = i // image_width
            y_max = y_min
            is_valid = True

            stack = [i]
            while stack:
                idx = stack.pop()
                if idx < 0 or idx >= total_pixels:
                    continue
                if image[idx] < cutoff or idx in checked:
                    continue

                if (idx % image_width == 0
                        or idx % image_width == image_width - 1
                        or idx // image_width == 0
                        or idx // image_width == image_height - 1):
                    is_valid = False

                checked.add(idx)
                star_indices.append(idx)

                if image[idx] > max_intensity:
                    max_intensity = image[idx]
                    guess = idx

                px = idx % image_width
                py = idx // image_width
                if px > x_max:
                    x_max = px
                elif px < x_min:
                    x_min = px
                if py > y_max:
                    y_max = py
                elif py < y_min:
                    y_min = py

                # Push in reverse order so RIGHT is processed first (LIFO),
                # matching C++ recursive right-first traversal
                stack.append(idx - image_width)  # UP (unconditional in C++)
                stack.append(idx + image_width)  # DOWN (unconditional in C++)
                if px != 0:
                    stack.append(idx - 1)        # LEFT
                if px != image_width - 1:
                    stack.append(idx + 1)        # RIGHT — last pushed = first popped

            x_diameter = (x_max - x_min) + 1

            if not is_valid:
                continue

            # Calculate FWHM
            count_fwhm = 0
            for idx in star_indices:
                if image[idx] > max_intensity / 2:
                    count_fwhm += 1
            fwhm = math.sqrt(count_fwhm)
            stddev = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            modified_stddev = 2.0 * stddev * stddev

            guess_x = float(guess % image_width)
            guess_y = float(guess // image_width)

            change = float('inf')
            stop = 0

            while change > self.MIN_CHANGE and stop < 100000:
                y_weighted_sum = 0.0
                x_weighted_sum = 0.0
                weighted_sum = 0.0
                stop += 1

                for idx in star_indices:
                    curr_x = float(idx % image_width)
                    curr_y = float(idx // image_width)
                    dx = curr_x - guess_x
                    dy = curr_y - guess_y
                    w = (max_intensity
                         * math.exp(-((dx * dx) / modified_stddev
                                      + (dy * dy) / modified_stddev)))

                    x_weighted_sum += w * curr_x * image[idx]
                    y_weighted_sum += w * curr_y * image[idx]
                    weighted_sum += w * image[idx]

                if weighted_sum > 0:
                    x_temp = x_weighted_sum / weighted_sum
                    y_temp = y_weighted_sum / weighted_sum
                    change = abs(guess_x - x_temp) + abs(guess_y - y_temp)
                    guess_x = x_temp
                    guess_y = y_temp
                else:
                    break

            result.append(Star(
                guess_x + 0.5, guess_y + 0.5,
                x_diameter / 2.0, y_diameter / 2.0,
                len(star_indices),
            ))

        return result
