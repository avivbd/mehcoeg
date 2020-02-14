classdef runTests < matlab.unittest.TestCase
    
    properties
            TestData
    end
        
    methods (Test)
                    
        function setupOnce(testCase)
            cd('../../../') % move back up to mehcoeg level
            testCase.TestData = load_data();
        end
                
        function testParams(testCase) 
            p = model_params();
            testCase.assertTrue(isa(p, 'struct'))
        end
        
    end
    
end 


