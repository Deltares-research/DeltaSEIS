import cv2
from pathlib import Path
import argparse

def resize_image(image, max_width):
    height, width = image.shape[:2]
    if width > max_width:
        scaling_factor = max_width / width
        new_width = max_width
        new_height = int(height * scaling_factor)
        return cv2.resize(image, (new_width, new_height))
    return image  # Return original image if no resizing is needed

def create_video_from_images(image_folder_path, video_name=None, max_width=1920, fps=30, cut_percentage=0):
    image_folder = Path(image_folder_path)
    images = sorted(image_folder.glob('*.jpg'))

    if video_name is None:
        video_name = image_folder.parent/"movie.avi"

    if not images:
        print(f"No images found in the specified folder: {image_folder_path}.")
        return

    total_images = len(images)
    cut_images = int(total_images * (cut_percentage / 100))

    # Adjust the images list to exclude the last cut_percentage
    images_to_process = images[:-cut_images] if cut_images > 0 else images

    if not images_to_process:
        print("No images left to process after cutting.")
        return

    first_image = cv2.imread(str(images_to_process[0]))
    first_image_resized = resize_image(first_image, max_width)
    height, width, layers = first_image_resized.shape

    print(f"Image dimensions after resizing: {width}x{height}")

    # Set video codec
    fourcc = cv2.VideoWriter_fourcc(*'XVID')  # or try 'MJPG'
    video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

    for image_path in images_to_process:
        img = cv2.imread(str(image_path))
        if img is None:
            print(f"Warning: Could not read image {image_path}. Skipping.")
            continue
        
        img_resized = resize_image(img, max_width)
        video.write(img_resized)

    video.release()
    cv2.destroyAllWindows()
    print(f"Video '{video_name}' created successfully!\n")

import argparse
from pathlib import Path

def create_video_from_images(image_folder_path, video_name=None, max_width=1920, fps=30, cut_percentage=0):
    if video_name is None:
        video_name = image_folder_path.parent / "movie.avi"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a video from JPG images.")
    parser.add_argument("image_folder", type=str, help="Path to the folder containing JPG images.")
    parser.add_argument("--video_name", type=str, help="Name of the output video file.")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second for the output video (default: 30).")
    parser.add_argument("--cut", type=float, default=0, help="Percentage of the video to cut from the end (default: 0).")
    
    args = parser.parse_args()
    
    # Convert image_folder to a Path object
    image_folder_path = Path(args.image_folder)
    
    # Call the function with the appropriate arguments
    create_video_from_images(image_folder_path, video_name=args.video_name, fps=args.fps, cut_percentage=args.cut)

