import numpy as np

def deploy_roi(center_lat_deg, center_lon_deg, radius_deg, segments):

    lats = []
    longs = []
    center_lat_ = 90 - center_lat_deg

    # Step size for generating points along the circle
    step = (2 * np.pi) / (segments - 1)
    theta = 0

    # Loop through each segment to compute the coordinates
    for i in range(segments):

        theta += step
        pointer_lat_ = np.arccos(np.cos(np.radians(radius_deg)) * np.cos(np.radians(center_lat_)) +
                                 np.sin(np.radians(radius_deg)) * np.sin(np.radians(center_lat_)) * np.cos(theta))
        lat = np.degrees(np.pi / 2 - pointer_lat_)

        # Longitude computation
        H = 1
        phi_mod_360 = np.mod(center_lat_, 360)
        if 180 <= phi_mod_360 < 360:
            H = -1

        delta_lon = np.arccos((np.cos(np.radians(radius_deg)) - np.cos(pointer_lat_) * np.cos(np.radians(center_lat_))) /
                              (H * np.sin(np.radians(center_lat_)) * np.sin(pointer_lat_))) - (np.pi / 2) * (H - 1)

        if theta <= np.pi:
            lon = center_lon_deg - np.degrees(delta_lon)
        else:
            lon = center_lon_deg + np.degrees(delta_lon)

        # Handling for poles
        if center_lat_deg == -90:
            lon = np.degrees(theta)
            lat = -90 + radius_deg
        elif center_lat_deg == 90:
            lon = np.degrees(theta)
            lat = 90 - radius_deg

        if theta == 2 * np.pi and (center_lat_deg + radius_deg) == 90:
            lat = 90
            lon = 0
        elif theta == 2 * np.pi and (center_lat_deg + radius_deg) == -90:
            lat = -90
            lon = 0

        # Correct longitude if needed
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360

        # Altitude is set to 0 for all points
        alt = 0

        # Ensure lat and lon are real numbers (sometimes precision errors cause a little imaginary part)
        lat = np.real(lat)
        lon = np.real(lon)

        # Append the calculated point to the list
        # roi_coordinates.append([lat, lon, alt])
        lats.append(lat)
        longs.append(lon)

    return lats,longs