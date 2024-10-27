import numpy as np

def compute_offset_points(points, distance):
    """
    Computes offset points at a specified distance perpendicular to the line segments
    defined by consecutive points in the input list.
    Args:
        points (list of tuple or list of list): A list of points where each point is 
            represented as a tuple or list of two coordinates (x, y).
        distance (float): The perpendicular distance at which to compute the offset points.
    Returns:
        list of numpy.ndarray: A list of offset points as numpy arrays. Each pair of 
        consecutive points in the input list will generate two offset points.
    """

    offset_points = []

    for i in range(len(points) - 1):
        p1 = np.array(points[i])
        p2 = np.array(points[i + 1])
        
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
        
        # Scale to the desired offset distance
        offset_vector = normal_vector * (distance / np.linalg.norm(normal_vector))
        
        # Calculate offset points
        offset_p1 = p1 + offset_vector
        offset_p2 = p2 + offset_vector
        
        offset_points.append(offset_p1)
        offset_points.append(offset_p2)
    
    return offset_points


