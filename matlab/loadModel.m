function loadModel
if ~libisloaded("model")
    loadlibrary('../Cpp/model.so','../Cpp/model.h')
end
