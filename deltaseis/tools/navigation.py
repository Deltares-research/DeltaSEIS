import numpy as np

def compute_offset_points(x, y, crossline_distance=2, inline_distance=0):
    """
    Offsets each point by a specified distance perpendicular and parallel to the vector between that point and the next point.
    Args:
        x (list of float): A list of x coordinates.
        y (list of float): A list of y coordinates.
        crossline_distance (float): The perpendicular distance at which to offset the points. Default is 2 meters.
        inline_distance (float): The parallel distance at which to offset the points. Default is 0 meters.
    Returns:
        tuple of lists: Two lists containing the x and y coordinates of the offset points.
    """

    x_offset = []
    y_offset = []

    for i in range(len(x) - 1):
        p1 = np.array([x[i], y[i]])
        p2 = np.array([x[i + 1], y[i + 1]])
        
        # Calculate direction vector
        direction = p2 - p1
        
        # Calculate the length of the direction vector
        length = np.linalg.norm(direction)
        
        if length == 0:
            continue  # Skip if points are the same
        
        # Normalize the direction vector
        direction_normalized = direction / length
        
        # Get the perpendicular vector (90 degrees rotation)
        normal_vector = np.array([-direction_normalized[1], direction_normalized[0]])
        
        # Scale to the desired crossline offset distance
        crossline_offset_vector = normal_vector * (crossline_distance / np.linalg.norm(normal_vector))
        
        # Scale to the desired inline offset distance
        inline_offset_vector = direction_normalized * inline_distance
        
        # Calculate offset point
        offset_point = p1 + crossline_offset_vector + inline_offset_vector
        
        x_offset.append(offset_point[0])
        y_offset.append(offset_point[1])
    
    # Handle the last point separately
    if len(x) > 1:
        p_last = np.array([x[-1], y[-1]])
        p_prev = np.array([x[-2], y[-2]])
        
        # Calculate direction vector
        direction = p_last - p_prev
        
        # Calculate the length of the direction vector
        length = np.linalg.norm(direction)
        
        if length != 0:
            # Normalize the direction vector
            direction_normalized = direction / length
            
            # Get the perpendicular vector (90 degrees rotation)
            normal_vector = np.array([-direction_normalized[1], direction_normalized[0]])
            
            # Scale to the desired crossline offset distance
            crossline_offset_vector = normal_vector * (crossline_distance / np.linalg.norm(normal_vector))
            
            # Scale to the desired inline offset distance
            inline_offset_vector = direction_normalized * inline_distance
            
            # Calculate offset point
            offset_point = p_last + crossline_offset_vector + inline_offset_vector
            
            x_offset.append(offset_point[0])
            y_offset.append(offset_point[1])
    
    return x_offset, y_offset


