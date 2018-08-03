clear all, close all, clc;

p = 2;     % base
Np = 4000;  % Nb projections
Ns = 10;   % Partition of the circle into Ns segments
Nz = 10000;   % Add a projection at a reference angle each Nz projections

Seq = Low_discrepency_squence(Np,Ns,p,Nz);
Seq_reshap = reshape(Seq,Ns,floor(Np/Ns))';



name = 'yflo';
filename = 'ctprofile.xml';
CtScannXML = parseXML(filename);

idx=1;

for ang = 1:floor(Np/Ns)

CtScannXML.Children(2).Children.Data =  [num2str(idx) '.' name];
idx = idx + 1;
CtScannXML.Children(4).Children(8).Children.Data = num2str(Seq_reshap(ang,1));

nbChildren = (size(CtScannXML.Children,2)-1)/2;

docNode = com.mathworks.xml.XMLUtils.createDocument(CtScannXML.Name);
docRootNode = docNode.getDocumentElement;

for i=1:size(CtScannXML.Attributes(),2)
    docRootNode.setAttribute(CtScannXML.Attributes(i).Name,CtScannXML.Attributes(i).Value);
end

for i=1:nbChildren
    thisElement = docNode.createElement(CtScannXML.Children(2*i).Name);
    
    CtScannXML.Children(2*i).Name
    
    nbGC = (size(CtScannXML.Children(2*i).Children,2)-1)/2;
    
    if nbGC<0
        docRootNode.appendChild(thisElement);
    elseif nbGC==0
        thisElement.appendChild(docNode.createTextNode(...
            CtScannXML.Children(2*i).Children.Data));
        docRootNode.appendChild(thisElement);
        
    else
    	for j = 1:nbGC
            thisChild = docNode.createElement(CtScannXML.Children(2*i).Children(2*j).Name);
            
            nbGGC = (size(CtScannXML.Children(2*i).Children(2*j).Children,2)-1)/2;
            
            if nbGGC==0
                thisChild.appendChild(docNode.createTextNode(...
                    CtScannXML.Children(2*i).Children(2*j).Children.Data));
                
            else
                for k = 1:nbGGC
                    thisGChild = docNode.createElement(CtScannXML.Children(2*i).Children(2*j).Children(2*k).Name);
                    thisGChild.appendChild(docNode.createTextNode(...
                        CtScannXML.Children(2*i).Children(2*j).Children(2*k).Children.Data));
                    thisChild.appendChild(thisGChild);
                end
            end
            thisElement.appendChild(thisChild);
            
            if size(CtScannXML.Children(2*i).Attributes(),2)
                for att=1:size(CtScannXML.Children(2*i).Attributes(),2)
                    thisElement.setAttribute(CtScannXML.Children(2*i).Attributes(att).Name,CtScannXML.Children(2*i).Attributes(att).Value);
                end
            end
            docRootNode.appendChild(thisElement);    
        end
    end
end

xmlFileName = strcat('XMLfiles/', name , sprintf('%04d',ang) ,'.ctprofile.xml');
xmlwrite(xmlFileName,docNode);
type(xmlFileName);
end
