#using Pkg
#Pkg.add("ImageView")
#Pkg.add("Images")
#Pkg.add("ImageIO")
#Pkg.add("ImageMagick")
#Pkg.add("TSVD")
using TSVD
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

basefilename="roundabout_snapshots"

img=load(abspath(basefilename*"_0001.png"));
# Visualize it
#imshow(img);
# Determine the image size
sz=size(img); szv=sz[1]*sz[2];
# Reshape the image to a vector:
R=float(red.(img));
G=float(green.(img));
B=float(blue.(img));
v=vcat(vec(R), vec(G), vec(B))
#
m=56  # Number of frames to load
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
U, s, V = TSVD.tsvd(A, 5);
#display(s);
A_rank_1 = U[:,1:2]*Diagonal(s[1:2])*V[:,1:2]';
#A_rank_1 = U[:,1]*s[1]*V[:,1]';

for i in 1:1
    v=A_rank_1[:,i]; # Third column of A = third frame
    vv=reshape(v,szv,3);
    R=reshape(vv[:,1],sz[1],sz[2]);
    G=reshape(vv[:,2],sz[1],sz[2]);
    B=reshape(vv[:,3],sz[1],sz[2]);
    rank1_img=RGB.(R,G,B);
    imshow(rank1_img);
    save("excercise_6c_rank2_$i.png",rank1_img);
    #save("excercise_6c_rank1_$i.png",rank1_img);
end
