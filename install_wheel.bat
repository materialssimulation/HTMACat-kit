pip install ase numpy==1.23.5 scikit-learn
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/FireWorks-2.0.3-py3-none-any.whl
pip install https://raw.githubusercontent.com/stanfordbshan/HTMACat-kit/master/requires_wheel/CatKit-0.5.4-py3-none-any.requires_wheel
for /f "tokens=1,* delims=:" %%A in ('curl -ks https://api.github.com/repos/stanfordbshan/HTMACat-kit/releases/latest ^| find "browser_download_url"') do (
    pip install -U %%B
)