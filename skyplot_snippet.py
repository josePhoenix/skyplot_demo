# -*- coding: utf-8 -*-
from __future__ import print_function
import glob
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import datetime
import os.path
import ephem

def equatorial_deg_to_ecliptic_deg(ra, dec):
    """Convert RA and declination as decimal degrees to
    ecliptic longitude and latitude in decimal degrees"""
    eq = ephem.Equatorial(ra * ephem.degree, dec * ephem.degree)
    ec = ephem.Ecliptic(eq)
    return ec.lon / ephem.degree, ec.lat / ephem.degree

def ecliptic_deg_to_equatorial_deg(lon, lat):
    """Convert ecliptic longitude and latitude as decimal degrees to
    equatorial RA and declination as decimal degrees"""
    ec = ephem.Ecliptic(lon * ephem.degree, lat * ephem.degree)
    eq = ephem.Equatorial(ec)
    return eq.ra / ephem.degree, eq.dec / ephem.degree

def equatorial_plane_rad():
    """Get ecliptic longitude and latitude coordinates in radians
    for the celestial equator (declination = 0 deg)"""
    equatorial_plane_ra = np.linspace(0, 360, 200)
    equatorial_plane_dec = np.zeros_like(equatorial_plane_ra)
    equatorial_plane_l, equatorial_plane_b = np.empty_like(equatorial_plane_ra), np.empty_like(equatorial_plane_dec)
    for idx, (ra, dec) in enumerate(zip(equatorial_plane_ra, equatorial_plane_dec)):
        l, b = equatorial_deg_to_ecliptic_deg(ra, dec)
        equatorial_plane_l[idx] = l
        equatorial_plane_b[idx] = b
    return np.deg2rad(equatorial_plane_l), np.deg2rad(equatorial_plane_b)

def galactic_plane_rad():
    """Get ecliptic longitude and latitude coordinates in radians
    for the galactic plane"""
    galactic_l = np.linspace(0, 2 * np.pi, 200)
    galactic_b = np.zeros_like(galactic_l)
    ecliptic_l, ecliptic_b = np.empty_like(galactic_l), np.empty_like(galactic_b)
    for idx, (g_l, g_b) in enumerate(zip(galactic_l, galactic_b)):
        gal = ephem.Galactic(g_l, g_b)
        ec = ephem.Ecliptic(gal)
        ecliptic_l[idx] = ec.lon
        ecliptic_b[idx] = ec.lat
    return ecliptic_l, ecliptic_b

def dest_from_start_and_bearing_rad(lon, lat, bearing, angular_dist):
    """
    Find a destination's latitude and longitude given a
    starting (longitude, latitude) pair, a bearing in
    [0, 2pi], and an angular distance. (All arguments in radians.)
    
    From the explanation at
    http://www.movable-type.co.uk/scripts/latlong.html#destPoint
    """
    lat2 = np.arcsin(np.sin(lat) * np.cos(angular_dist) + np.cos(lat) * np.sin(angular_dist) * np.cos(bearing))
    lon2 = lon + np.arctan2(np.sin(bearing) * np.sin(angular_dist) * np.cos(lat), np.cos(angular_dist) - np.sin(lat) * np.sin(lat2))
    return lon2, lat2

def small_circle_rad(lon, lat, angular_radius, npoints=100):
    """Compute longitude and latitude coordinates for `npoints`
    points in a circle around longitude `lon` and latitude `lat` in radians"""
    lons = np.ones(npoints) * lon
    lats = np.ones(npoints) * lat
    angular_radii = np.ones(npoints) * angular_radius
    bearings = np.linspace(0, 2 * np.pi, npoints)
    lonpts, latpts = dest_from_start_and_bearing_rad(lons, lats, bearings, angular_radii)
    return lonpts, latpts

def _plot_deg(func, lon, lat, **kwargs):
    """Take longitude degrees [0, 360] to radians [-pi, pi], latitude degrees to radians, and plot"""
    if not np.isscalar(lon):
        lon = lon.copy()
        lon %= 360.0
        lon[lon > 180] -= 360
        lon_rad, lat_rad = np.deg2rad(lon), np.deg2rad(lat)
    else:
        lon %= 360.0
        lon = lon - 360 if lon > 180 else lon
        lon_rad, lat_rad = np.deg2rad(lon), np.deg2rad(lat)
    func(lon_rad, lat_rad, **kwargs)

def _plot_rad(func, lon, lat, **kwargs):
    """Wrap radians at longitude pi and plot"""
    if not np.isscalar(lon):
        lon = lon.copy()
        lon %= 2 * np.pi
        lon[lon > np.pi] -= 2 * np.pi
        lon, lat = nanify(lon, lat)
    else:
        lon %= 2 * np.pi
        lon = lon - 2 * np.pi if lon > np.pi else lon
    func(lon, lat, **kwargs)

def plot_deg(ax, lon, lat, **kwargs):
    _plot_deg(ax.plot, lon, lat, **kwargs)
def scatter_deg(ax, lon, lat, **kwargs):
    _plot_deg(ax.scatter, lon, lat, **kwargs)

def plot_rad(ax, lon, lat, **kwargs):
    _plot_rad(ax.plot, lon, lat, **kwargs)
def scatter_rad(ax, lon, lat, **kwargs):
    _plot_rad(ax.scatter, lon, lat, **kwargs)

def nanify(lon, lat):
    """Take longitude and latitude in radians, replace
    certain `lon` values where coordinates wrap around
    with `np.nan` to create discontinuities
    (and prevent weird wrapping behavior)"""
    wjump = np.where(((lon[:-1] >  np.pi/2) & (lon[1:] < -np.pi/2)) | 
                     ((lon[:-1] < -np.pi/2) & (lon[1:] >  np.pi/2)))
    lat = lat.copy()
    lat[wjump] = np.nan

    return lon, lat

sun_angle_from, sun_angle_to = 85, 85 + 50  # can point up to 5 deg towards sun, up to 50 deg away from sun

def check_separation(point_a, point_b, radius_degrees):
    """
    ((long, lat) in deg, (long, lat) in deg, radius in deg) -> True or False

    If the central angle between the two points, as computed by the haversine
    formula, is less than or equal to the radius in degrees, True is returned.

    https://en.wikipedia.org/wiki/Haversine_formula
    """
    import math

    def haversine(theta):
        return math.sin(theta / 2.0) ** 2

    hav_d_over_r = haversine(phi2 - phi1) + \
        math.cos(phi1) * math.cos(phi2) * haversine(lambda2 - lambda1)

    central_angle_rad = 2 * math.asin(math.sqrt(hav_d_over_r))
    central_angle_deg = central_angle_rad * 180.0 / math.pi

    if central_angle_deg <= radius_degrees:
        return True
    else:
        return False

def angular_separation_rad(lambda1, phi1, lambda2, phi2):
    """
    Compute the angular separation from (lambda1, phi1)
    in radians to (lambda2, phi2) in radians.
    
    lambda1, phi1 : arrays or scalars
    lambda2, phi2 : scalars
    """
    def haversine(theta):
        return np.sin(theta / 2.0) ** 2
    hav_d_over_r = haversine(phi2 - phi1) + \
        np.cos(phi1) * np.cos(phi2) * haversine(lambda2 - lambda1)
    central_angle_rad = 2 * np.arcsin(np.sqrt(hav_d_over_r))
    return central_angle_rad
    
def field_of_regard_filter(el, eb, sun):
    lambdas, phis = np.deg2rad(el), np.deg2rad(eb)
    separations = angular_separation_rad(lambdas, phis, sun.lon, sun.lat)
    sun_angle_from_rad, sun_angle_to_rad = np.deg2rad(sun_angle_from), np.deg2rad(sun_angle_to)
    rows_with_separation_in_range = (separations > sun_angle_from_rad) & (separations < sun_angle_to_rad)

    return rows_with_separation_in_range

def plot_sky(sky_ax, el, eb, for_date=None):
    """
    Plot the sky and a catalog in ecliptic coordinates
    
    Parameters
    ----------
    sky_ax : Axes instance
        Axes (initialized using the mollweide projection)
        on which to plot
    el, eb : arrays or sequences
        Ecliptic longitude (`el`) and latitude (`eb`) coordinates for
        the catalog in degrees
    for_date : datetime.date instance, optional
        Use a given date to compute the sun position. Default is today.
    """
    if for_date is None:
        for_date = datetime.date.today()
    el, eb = np.asarray(el), np.asarray(eb)
    sky_ax.set_title('Ecliptic Coordinates')
    sky_ax.grid()
    
    # locate sun
    sun = ephem.Sun()
    sun.compute(for_date)
    sun_ec = ephem.Ecliptic(sun)

    
    subset_mask = field_of_regard_filter(el, eb, sun_ec)
    # plot all points
    scatter_deg(sky_ax, el, eb, marker='.', alpha=0.5, color='blue')
    # plot points ~within~ FoR
    scatter_deg(sky_ax, el[subset_mask], eb[subset_mask], marker='.', color='red', alpha=0.5)

    # plot sun
    plot_deg(sky_ax, sun_ec.lon / ephem.degree, sun_ec.lat / ephem.degree, markersize=15, marker='o', color='yellow')
    # ... and anti-sun
    plot_deg(sky_ax, sun_ec.lon / ephem.degree + 180, 0, markersize=15, marker='o', color='red')

    # add equatorial plane
    equatorial_plane_l, equatorial_plane_b = equatorial_plane_rad()
    plot_rad(sky_ax, equatorial_plane_l, equatorial_plane_b, label='Equator', lw=2, alpha=0.5)

    # add galactic plane
    gal_l, gal_b = galactic_plane_rad()
    plot_rad(sky_ax, gal_l, gal_b, lw=2, label='Galactic', alpha=0.5)
    
    # plot field of regard
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_from))
    plot_rad(sky_ax, cl, cb, lw=4, label=u"{}º to {}º from sun".format(sun_angle_from, sun_angle_to), color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad(sun_angle_to))
    plot_rad(sky_ax, cl, cb, lw=4, color='orange')
    cl, cb = small_circle_rad(sun_ec.lon, sun_ec.lat, np.deg2rad((sun_angle_to + sun_angle_from) / 2.0))
    plot_rad(sky_ax, cl, cb, lw=4, ls='--', color='orange')
    # sky_ax.legend()
