#using Pkg
#Pkg.add("ImageView")
#Pkg.add("Images")
#Pkg.add("ImageIO")
#Pkg.add("ImageMagick")
#Pkg.add("FileIO")
#Pkg.add("TestImages")
using ImageView
using Images
using ImageIO
using ImageMagick
using Printf
using FileIO
using LinearAlgebra
using Arpack
using TestImages
using Plots

basefilename="india_driving_frame"

img=load(abspath(basefilename*"_0001.png"));
# Visualize it
# Determine the image size
sz=size(img); szv=sz[1]*sz[2];
# Reshape the image to a vector:
R=float(red.(img));
G=float(green.(img));
B=float(blue.(img));
v=vcat(vec(R), vec(G), vec(B))
#
m=3  # Number of frames to load
A=zeros(size(v,1),m) # Matrix to store all the frames in
for k=1:m
    fname=@sprintf("%s_%04d.png",basefilename,k);
    img=load(abspath(fname));
    #println(fname)
    R=float(red.(img));
    G=float(green.(img));
    B=float(blue.(img));
    v=vcat(vec(R), vec(G), vec(B));
    A[:,k]=v;
end
U,σ,V = svd(A);
display(σ);
#A_rank_1 = U[:,1:3]*Diagonal(σ[1:3])*V[:,1:3]';
#A_rank_1 = U[:,1]*σ[1]*V[:,1]';

#for i in 1:1
#    v=A_rank_1[:,i]; # Third column of A = third frame
#    vv=reshape(v,szv,3);
#    R=reshape(vv[:,1],sz[1],sz[2]);
#    G=reshape(vv[:,2],sz[1],sz[2]);
#    B=reshape(vv[:,3],sz[1],sz[2]);
#    rank1_img=RGB.(R,G,B);
#    imshow(rank1_img);
#    save("excercise_6b_rank1_$i.png",rank1_img);
#end
