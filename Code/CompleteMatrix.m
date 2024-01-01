function completed_matrix = CompleteMatrix(matrix)
        mask = isnan(matrix);
        matrix(isnan(matrix)) = 0;
        

        [completed_matrix,~] = MC_Nuclear_IALM(matrix,~mask);
        
  
end