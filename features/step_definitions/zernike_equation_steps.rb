require 'uri'
require 'cgi'
require File.expand_path(File.join(File.dirname(__FILE__), "..", "support", "paths"))
#require File.expand_path(File.join(File.dirname(__FILE__), "..", "support", "selectors"))

module WithinHelpers
  def with_scope(locator)
    locator ? within(*selector_for(locator)) { yield } : yield
  end
end
World(WithinHelpers)


When /^(?:|I )go to "(.+)"$/  do |page_name|
    visit('/zernikes')
    # visit path_to(page_name)
end

Then /^I click "(.+)"$/ do |link|
    click_link(link)
end

When /^(?:|I )fill in "([^"]*)" with "([^"]*)"$/ do |field, value|
  fill_in(field, :with => value)
end

When /^(?:|I )fill in "([^"]*)" for "([^"]*)"$/ do |value, field|
  fill_in(field, :with => value)
end


When /^(?:|I )press "(.*)"$/ do |button|
  click_button(button)
end

Then /^I should see "(.*)"$/ do |text|
    if page.respond_to? :should
        page.should have_content(text)
    else
        assert page.has_content?(text)
    end
end


                                                                                                         
