% testing routine for writing triple prod

function [] = test_write_triple_prod(M_max)

  error_Dxi = [];
  error_D2xi = [];
  
  for m = 2 : M_max
   filename_Dxi = strcat('triple_prod/triple_prod_Dxi_He',num2str(m),'.txt');
   filename_D2xi = strcat('triple_prod/triple_prod_D2xi_He',num2str(m),'.txt');
   
   read_Dxi = dlmread(filename_Dxi,'\t');
   read_Dxi = sparse(read_Dxi(:,1),read_Dxi(:,2),read_Dxi(:,3));

   read_D2xi = dlmread(filename_D2xi,'\t');
   read_D2xi = sparse(read_D2xi(:,1),read_D2xi(:,2),read_D2xi(:,3));
   
   manual_mat = 0 * read_Dxi;
   manual_mat2 = 0 * read_D2xi;
   
   for i = 1 : size(manual_mat,1)
       for j = 1 : size(manual_mat,2)
           manual_mat(i,j) = He_triple_product_Dxi(m,i-1,j-1);
           manual_mat2(i,j) = He_triple_product_D2xi(m,i-1,j-1);
       end
   end
   
   error_Dxi = [norm(full(read_Dxi)-manual_mat), error_Dxi];
   error_D2xi = [norm(full(read_D2xi)-manual_mat2), error_D2xi];
  end
  
  disp('error Dxi');
  norm(error_Dxi)
  
  disp('error D2xi');
  norm(error_D2xi)
end