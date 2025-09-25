import numpy as np
import matplotlib.pyplot as plt

from matplotlib.image import imread
from PIL import Image, ImageDraw, ImageFont
import os


# os.makedirs('results')

PATH = os.path.join(os.getcwd(),'data')

def make_folder(folder_name:str,path:str):

    if folder_name not in os.listdir(path):
        os.makedirs(os.path.join(path,folder_name))
    
    file_path = os.path.join(path,folder_name)
    
    return file_path
    

def make_plots(filename:str):
    file_path = f'{PATH}\\{filename}.dat'
    
    data = np.genfromtxt(file_path,skip_header=1,delimiter=' ')
    label = open(file_path).readline().strip().split(' ')
    
    x = data[:,0]
    y = data[:,1]

    fig,ax = plt.subplots()
    fig.set_size_inches((10,6))
    
    ax.set_xlabel(f'{label[0]}')
    ax.set_ylabel(f'{label[1]}')
    ax.plot(x,y)
    
    output_path = make_folder('results',PATH)
    
    fig.savefig(f'{output_path}\\{filename}.png')
    
    

# def make_gif(foldername:str,gif_name:str):
    
#     folder_path = f'{PATH}\\{foldername}'
#     output_path = make_folder('results',PATH)
    
    
#     imgs = [i for i in os.listdir(folder_path) if i.endswith('.ppm')]
#     all_images = []

#     for file in imgs:
#         image_path = os.path.join(folder_path,file)
        
#         path = 'density'
#         font = ImageFont.truetype(font="arial.ttf", size=20)

#         for img in os.listdir(path)[:1]:

#             image_path = os.path.join(path,file)
#             image = Image.open(image_path)
#             draw = ImageDraw.Draw(image)
            
#             draw.text((5,5),text=f"Iteration :\n{img.strip('.ppm').strip(path)}",font=font)
#             all_images.append(Image.open(image_path))
        
#     all_images[0].save(f'{output_path}\\{gif_name}.gif', save_all = True ,append_images = all_images[1:])
    
def make_gif(foldername:str,gif_name:str):
    print(f"creating the gif for {gif_name}")
    folder_path = f'{PATH}\\{foldername}'
    output_path = make_folder('results',PATH)


    imgs = [i for i in os.listdir(folder_path) if i.endswith('.ppm')]
    all_images = []
    font = ImageFont.truetype(font="arial.ttf", size=20)

    for file in imgs:
        image_path = os.path.join(folder_path,file)
        image = Image.open(image_path)
        draw = ImageDraw.Draw(image)
        draw.text((5,5),text=f"Iteration :\n{file.strip('.ppm').strip(foldername)}",font=font)

        all_images.append(image)
        
    all_images[0].save(f'{output_path}\\{gif_name}.gif', save_all = True ,append_images = all_images[1:])
    print(f"{gif_name}.gif is saved at {output_path}\\{gif_name}.gif ")



def plots_with_fields(filename:str,foldername:str,output_name:str):
     
    file_path = f'{PATH}\\{filename}.dat'
    folder_path = f'{PATH}\\{foldername}'
    
    
    label = open(file_path).readline().strip().split(' ')
    data = np.genfromtxt(file_path,skip_header=1,delimiter=' ')
    
    output_path = make_folder('results',PATH)
    
    imgs = [i for i in os.listdir(folder_path) if i.endswith('.ppm')][-1]
    
    image_path = os.path.join(folder_path,imgs)
    shape_fr = (200, 60)

    ux = data[:, 2].reshape(shape_fr)
    uy = data[:, 3].reshape(shape_fr)
        # Flattened coordinates
    x_flat = data[:, 0]
    y_flat = data[:, 1]

    # Get unique x and y values
    x_unique = np.unique(x_flat)
    y_unique = np.unique(y_flat)

    # Create meshgrid
    X, Y = np.meshgrid(x_unique, y_unique)

    # Reshape ux and uy to match meshgrid
    UX = ux.reshape(X.shape)
    UY = uy.reshape(Y.shape)


    fig,ax = plt.subplots()
    fig.set_size_inches((10,6))
    
    ax.set_xlabel(f'{label[0]}')
    ax.set_ylabel(f'{label[1]}')
    

    ax.streamplot(X, Y, UX, UY, color='k', cmap="plasma")
    ax.imshow(imread(image_path))
    
    fig.savefig(f'{output_path}\\{output_name}.png')
    
    pass

if __name__=='__main__':


    # make_plots('data')
    # make_gif(foldername='rho_heavy',gif_name='rho_heavy')
    make_gif(foldername='rho_light',gif_name='rho_light')

    # plots_with_fields(filename='data',foldername='vorticity',output_name='vorticity')
    # make_gif(foldername='vorticity',gif_name='vorticity')
    # plots_with_fields(filename='data',foldername='velocity',output_name='velocity')
    
    
    pass
