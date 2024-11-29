import asyncio
import aiohttp
import time
import numpy as np
import pandas as pd
from os.path import join as opj
from pathlib import Path
from asyncio.exceptions import TimeoutError, CancelledError


def convert_to_number(rn: str) -> int:
    return int(rn.replace("-", ""))


async def fetch_casrn_for_smiles(
    smiles: str,
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    idx: int,
    total: int,
) -> str:
    """Fetch CASRN for a given SMILES string from PubChem API asynchronously, with progress printing."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/xrefs/RN/JSON"

    async with semaphore:  # Ensure the number of simultaneous requests is limited
        try:
            async with session.get(
                url, timeout=20
            ) as response:  # Increased timeout to 20 seconds
                if response.status == 200:
                    json_result = await response.json()
                    casrns = json_result["InformationList"]["Information"][0]["RN"]
                    casrn = min(casrns, key=convert_to_number)
                    print(
                        f"{smiles} -> {casrn} ({idx + 1}/{total})"
                    )  # Progress printing
                    return str(convert_to_number(casrn))
        except aiohttp.ClientError as e:
            print(
                f"Error fetching CASRN for SMILES: {smiles}, aiohttp error: {e} ({idx + 1}/{total})"
            )
            return None
        except TimeoutError as e:
            print(f"TimeoutError: {smiles} timed out after waiting ({idx + 1}/{total})")
            return None  # Handle timeout specifically
        except CancelledError as e:
            print(f"CancelledError: {smiles} request was cancelled ({idx + 1}/{total})")
            return None  # Handle cancelled task
        except Exception as e:
            print(f"Unexpected error: {e} ({idx + 1}/{total})")
            return None  # Catch any other unexpected errors


async def process_smiles(smiles_list, max_concurrent_requests=5):
    """Process all SMILES using async requests with a limit on concurrent requests."""
    casrn_list = []
    failed_smiles = []
    semaphore = asyncio.Semaphore(
        max_concurrent_requests
    )  # Limit the number of concurrent requests
    total_smiles = len(smiles_list)  # Total SMILES for progress

    async with aiohttp.ClientSession() as session:
        tasks = [
            fetch_casrn_for_smiles(smiles, session, semaphore, idx, total_smiles)
            for idx, smiles in enumerate(smiles_list)
        ]
        results = await asyncio.gather(
            *tasks, return_exceptions=True
        )  # Allow exceptions in the result

    for i, (smiles, result) in enumerate(zip(smiles_list, results)):
        if result is not None:
            casrn_list.append(result)
        else:
            failed_smiles.append(smiles)
            casrn_list.append("")  # Append empty string for failed cases

    return casrn_list, failed_smiles


async def retry_failed_smiles(failed_smiles, max_retries=5, max_concurrent_requests=5):
    """Retry the failed SMILES requests up to max_retries."""
    semaphore = asyncio.Semaphore(max_concurrent_requests)
    retry_count = 0
    total_failed = len(failed_smiles)  # Track the total number of retries

    while failed_smiles and retry_count < max_retries:
        retry_count += 1
        print(f"Retrying failed SMILES... Attempt {retry_count}/{max_retries}")

        casrn_list = []
        new_failed_smiles = []
        async with aiohttp.ClientSession() as session:
            tasks = [
                fetch_casrn_for_smiles(smiles, session, semaphore, idx, total_failed)
                for idx, smiles in enumerate(failed_smiles)
            ]
            results = await asyncio.gather(
                *tasks, return_exceptions=True
            )  # Allow exceptions in the result

        for smiles, result in zip(failed_smiles, results):
            if result is not None:
                casrn_list.append(result)
                print(f"Retry success: {smiles} -> {result}")
            else:
                new_failed_smiles.append(smiles)
                casrn_list.append("")  # Append empty string for failed cases

        failed_smiles = new_failed_smiles

    return casrn_list, failed_smiles


def main():
    # Load the data
    df = pd.read_csv(opj(Path(__file__).parent.parent, "data", "DDB_Binary_VLE05.csv"))

    # Combine and deduplicate the SMILES lists
    smiles_list = list(set(df["cmp1_SMILES"].to_list() + df["cmp2_SMILES"].to_list()))

    # Run the async event loop to process the SMILES
    casrn_list, failed_smiles = asyncio.run(process_smiles(smiles_list))

    # Retry failed requests
    if failed_smiles:
        retry_casrn_list, final_failed_smiles = asyncio.run(
            retry_failed_smiles(failed_smiles)
        )
        # Replace failed entries with retry results
        for i, result in enumerate(retry_casrn_list):
            if result:  # Only update if there's a result
                index = smiles_list.index(failed_smiles[i])
                casrn_list[index] = result

    # Save the result to CSV and Excel with the "_refacv2" suffix
    output_df = pd.DataFrame(data={"CASRN": casrn_list, "SMILES": smiles_list})
    output_df.to_csv("CASRN_SMILES_refacv2.csv", index=False)
    output_df.to_excel("CASRN_SMILES_refacv2.xlsx", index=False)


if __name__ == "__main__":
    main()
