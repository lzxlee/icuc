Feature: inputting parameters into the Zernike equation
 
  As an optometry student
  I want to be able input parameters into Zernike equation
  So I can generate the Zernike equation

Scenario: manually inputting Zernike equation parameters
  When I go to "the home page"
  And I click "enter manually"
  And I press "Compute!"
  Then I should see "Zernike equation successfully computed!"

# Scenario: uploading Zernike equation parameters from a corrupted file
#   # When I go to "home"
#   When I go to "the home page"
#   And I click "upload text file"
#   Then I should see "Unable to upload file"