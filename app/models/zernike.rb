class Zernike < ActiveRecord::Base
  def self.parameters
    par = []
    (1..65).each do |i| 
        par << [i, "parameter #{i}"]
    end
    return par
  end
  
  
  def create
    
  end
end 